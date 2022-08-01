import os
from re import S
import utils
import plottools as pt
import pandas as pd
import numpy as np
from dateutil.relativedelta import relativedelta as relativedelta


GAS_ORDER=["co","so2","o3","no2"]
SENSOR_ORDER=[
            "c1_0#ec_na#0",
            "c1_5#ec_na#0",
            "c1_4#ec_na#0",
            "c1_3#ec_na#0"
            ]
K_=[ 
    [
        1,
        1.1100000143051147,
        0,
        3.630000114440918
    ],
    [
        0.12999999523162842,
        1,
        0,
        4.367000102996826
    ],
    [
        0,
        0,
        1,
        0
    ],
    [
        0,
        -0.15399999916553497,
        0,
        1
    ],
]

def corrected_ec_regression(point_surv_i:np.ndarray,point_surv_j:np.ndarray,point_stat_i:np.ndarray,point_stat_j:np.ndarray):
    k=np.asarray(K_)
    s_i=np.diag(point_surv_i)
    s_j=np.diag(point_surv_j)
    g_i=point_stat_i.reshape(-1,1)
    g_j=point_stat_j.reshape(-1,1)
    # k dot g =   A' dot s'  +B  =  [ s'1*A'1 , s'2*A'2,..., s'n*A'n  ]T +B = s dot A +B
    # k is n*n
    # g is n*1 (濃度)
    # s' is n*1  (電流)
    # s is n*n diagonal  (電流)
    # A' is n*n diagonal
    # A is n*1
    # B is n*1
    # ti:  k dot g_i =   s_i dot A +B
    # tj:  k dot g_j =   s_j dot A +B
    r_ti= np.dot(k,g_i)
    r_tj= np.dot(k,g_j)
    dr= r_ti-r_tj
    dks=s_i-s_j
    # A= dks.Inver * dr
    # B = r_ti - ks_i* A 

    dinv = np.linalg.inv(dks)
    A= np.dot(dinv,dr)
    B= r_ti - np.dot(s_i,A)
    return A.reshape(-1),B.reshape(-1)

def correction_all_points_processing(surv_vender:utils.RefSurveiesVender,stat_vender:utils.RefStationsVender,regression_date_list):
    #此為測站按時間的之目標小時aqi, 已為平均
    stat_regression_point_list= stat_vender.getDataByDateList(regression_date_list)
    #此為機器按時間的之一小時內所有點, 應取平均或某種取點法
    surv_regression_point_list= surv_vender.getDataByDateList(regression_date_list) 

    result=[]
    for current_date in regression_date_list:
        pass_points=utils.getIndex_of_current_date_and_before(current_date,regression_date_list)
        pass_points.reverse()
        if len(pass_points) <=1:
            continue
        #先檢查當前點
        point_index_i=pass_points[0]
        try:
            #此範例取最後一點
            surv_selected=surv_regression_point_list[point_index_i].iloc[-1]
            point_surv_i = surv_vender.get_value_array_of_na(SENSOR_ORDER,surv_selected)
        except utils.MissingValueException:
            continue #在這個範例 當前點 有缺漏值就不做校正
        
        try:
            point_stat_i = stat_vender.get_value_array_of_ppb(GAS_ORDER,stat_regression_point_list[point_index_i])
        except utils.MissingValueException:
            continue #在這個範例 當前點 有缺漏值就不做校正

        #檢查上一個點
        point_surv_j=None
        point_stat_j=None
        for index in pass_points[1:3]: #在這個範例至多只檢查就近兩點
            try:
                #此範例取最後一點
                surv_selected=surv_regression_point_list[index].iloc[-1]
                point_surv_j=surv_vender.get_value_array_of_na(SENSOR_ORDER,surv_selected)
                if surv_vender.check_diff_zreo(point_surv_i,point_surv_j):
                    point_surv_j=None
                    point_stat_j=None
                    continue #有零值就下一個
                    
                     
            except utils.MissingValueException:
                point_surv_j=None
                point_stat_j=None
                continue #有缺漏值就下一個
            try:
                point_stat_j=stat_vender.get_value_array_of_ppb(GAS_ORDER,stat_regression_point_list[index])
            except utils.MissingValueException:
                point_surv_j=None
                point_stat_j=None
                continue #有缺漏值就下一個
            
        if point_surv_j is None or point_stat_j is None:
            continue #檢查點皆缺漏則此次不校正
        
        a,b=corrected_ec_regression(point_surv_i,point_surv_j,point_stat_i,point_stat_j)
        result.append({"a":a,"b":b,"date":current_date})
        
    return result

    
def ec_regression(na_array:np.ndarray,a:np.ndarray,b:np.ndarray):
    # k dot g =   A' dot s'  +B  =  [ s'1*A'1 , s'2*A'2,..., s'n*A'n  ]T +B = s dot A +B
    # k is n*n
    # g is n*1 (濃度)
    # s' is n*1  (電流)
    # s is n*n diagonal  (電流)
    # A' is n*n diagonal
    # A is n*1
    # B is n*1
    # ti:  k dot g_i =   s_i dot A +B
    # tj:  k dot g_j =   s_j dot A +B
    k=np.asarray(K_)
    k_in=np.linalg.inv( k)
    # A dot I
    adi= np.dot(np.diag(a),na_array)
    g= np.dot(k_in ,adi) + np.dot(k_in ,b)
    return g

def ec_regression_all_points_processing(sync_points,surv_vender):
    current_sync_point_index=0
    surv_df=surv_vender.getData()
    def fn(row):
        nonlocal current_sync_point_index
        between_start=utils.any_datetime_2_utc_timestamp(sync_points[current_sync_point_index]["date"])
        if current_sync_point_index+1 < len(sync_points):
            between_end=utils.any_datetime_2_utc_timestamp(sync_points[current_sync_point_index+1]["date"])
        else:
            between_end=np.inf

        # 計算開始第二次校正前的部分
        if current_sync_point_index == 0 and row["m_time"] < between_start: 
            d={}
            for i,c in enumerate( GAS_ORDER):
                d[f"pred_{c}#ppb"]=np.NaN
            return d

        # 計算開始第二次校正後~最後一次校正
        if row["m_time"] >= between_start and row["m_time"] < between_end:
            sensor_na=surv_vender.get_value_array_of_na(SENSOR_ORDER,row)
            a=sync_points[current_sync_point_index]["a"]
            b=sync_points[current_sync_point_index]["b"]
            g_array=ec_regression(sensor_na,a,b)
            d={}
            for i,c in enumerate( GAS_ORDER):
                d[f"pred_{c}#ppb"]=g_array[i]
            return d

        # 該推進使用下一個點
        if row["m_time"] >= between_end:
            current_sync_point_index+=1
            return fn(row)
       

    df2=pd.DataFrame.from_records(surv_df.apply(fn,axis=1))
    return pd.concat([surv_df,df2], join='outer',axis=1)



def plot_full(surv_df:pd.DataFrame,stat_df:pd.DataFrame):
    pred_cols=[ f"pred_{c}#ppb" for c in GAS_ORDER]
    stat_cols=[ f"{c}#ppb" for c in GAS_ORDER]
    sensor_na_cols=SENSOR_ORDER
    import matplotlib.pyplot as plt  

    
    stat_df.plot(x='monitordate',y=stat_cols)
    plt.show()
    surv_df.plot(x='measure_time',y=pred_cols)
    plt.show()




if __name__ == '__main__':
    key_path = os.path.join('../tabairbox-dev-dydb-reader.csv')
    #機器本地庫存
    ref_surveies_pool_path='deviceData'
    #測站本地庫存
    ref_stations_pool_path='3rdData'
    station_name='aqx_p_190'
    device_uuid='69bdc804-30d0-43c7-a816-19f60e137952' #ETACOMA19000002
    #設定掃描時間區間與現在時間
    now=utils.datetime_utc_now()
    m_earlier=utils.any_isoformat_2_utc_datetime("2022-06-01 00:10:00Z")  #Z is UTC
    m_later=utils.any_isoformat_2_utc_datetime("2022-06-04 23:59:59Z") #Z is UTC
    #產生校正計算點時間清單 # https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.date_range.html
    regression_date_list=pd.date_range(m_earlier+relativedelta(hours=3),m_later,freq='1H').to_pydatetime() 

    #實例化測站庫存
    stat_storage=utils.RefStationsStorage(ref_stations_pool_path)
    #實例化測站提取工具
    stat_vender=utils.RefStationsVender(stat_storage,station_name,now,m_earlier,m_later)
    #實例化機器庫存
    surv_storage=utils.RefSurveiesStorage(key_path,ref_surveies_pool_path)
    #實例化機器提取工具
    surv_vender=utils.RefSurveiesVender(surv_storage,device_uuid,now,m_earlier,m_later)
    

    #對校正計算點時間清單進行計算,  得到各時間點a,b
    sync_points=correction_all_points_processing(surv_vender,stat_vender,regression_date_list)
    
    #進行濃度預測
    full_df=ec_regression_all_points_processing(sync_points,surv_vender)


    #畫圖

    stat_df = stat_vender.getData()

    station = pt.PlotStation(stat_df)
    station.plot_single()
    
    device_concentration = pt.PlotDeviceConcentration(full_df)
    device_concentration.plot_single()
    # device_concentration.plot_temp()
    
    device_current = pt.PlotDeviceCurrent(full_df)
    device_current.plot_single()


    





