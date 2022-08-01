import os,json,sys,csv,zipfile,re
import numpy as np
import pandas as pd
from pynamodb.models import Model as dynamoModel
from pynamodb.attributes import (
    UnicodeAttribute, NumberAttribute, UnicodeSetAttribute, BooleanAttribute,MapAttribute,JSONAttribute,ListAttribute
)
from dateutil.relativedelta import relativedelta as relativedelta
from dateutil.parser import parse
import datetime
import time,string,random
from typing import List, TypeVar,Callable,Union,Type,Tuple,Any,Dict
from typing import Callable
from typing import Tuple,Union
import pytz
from dateutil import tz


class MissingValueException(Exception):
    pass

def any_datetime_2_utc_timestamp(timeN:datetime.datetime)->float :
    if timeN.tzinfo is None:
        strptime=time.mktime(timeN.replace(microsecond=0).timetuple())
        return strptime
    else:
        timeN=timeN.astimezone(pytz.utc).replace(microsecond=0)
        return timeN.timestamp()


def utc_timestamp_2_utc_datetime(timeN:float)->datetime.datetime:
    return datetime.datetime.utcfromtimestamp( timeN).replace(tzinfo=pytz.utc)

def any_datetime_2_utc_isoformat(timeN:datetime.datetime)->str:
    return timeN.astimezone(pytz.utc).isoformat(timespec='milliseconds').replace('+00:00', 'Z')

def any_isoformat_2_utc_datetime(timeN:str)->datetime.datetime:
    #return parse_datetime.astimezone(pytz.utc)
    return parse(timeN).astimezone(pytz.utc)

def china_isoformat_2_utc_datetime(timeN:str)->datetime.datetime:
    #return parse_datetime.astimezone(pytz.utc)
    return parse(timeN).replace(tzinfo=datetime.timezone(datetime.timedelta(hours=8))).astimezone(pytz.utc)

def utc_timestamp_2_utc_isoformat(timeN:float)->str:
    return any_datetime_2_utc_isoformat(utc_timestamp_2_utc_datetime(timeN))

def any_isoformat_2_utc_timestamp(timeN:str)->float:
    return any_datetime_2_utc_timestamp(any_isoformat_2_utc_datetime(timeN))

def datetime_utc_now()->datetime.datetime:
    return datetime.datetime.utcnow().replace(tzinfo=pytz.utc).replace(microsecond=0)
def timestamp_utc_now()->float:
    #return round(time.time() * 1000)/1000
    return any_datetime_2_utc_timestamp(datetime_utc_now())


def get_sesnor_name_type(name):
    n=name.split('#')
    if len(n) <=2:
        return None
    return '{}#{}'.format(n[0], n[1])
def get_sesnor_name(name):
    n=name.split('#')
    if len(n) <=2:
        return None
    return n[0]

def get_sesnor_unit(name):
    n=name.split('#')
    if len(n) <=2:
        return None
    return n[1]

MOLECULAR_MASS= {
    'so2':64.06,  
    'no2':46.01,
    'co':28.01,
    'o3':48,
    'ethanol':46.07,
    'h2s':34.08,
    'no':30.01,
    'nh3':17.03,
    'co2':44.01,
    'cl2':70.90,
    'ch2o':30.03,
    'h2':2.02,
    'ch4':16.04,
    'cscl4':165.82,
    'voc':78.95
    }

def mcToVmd(input_unit_const,output_unit_const):
    #公式：重量體積密度 (mg/m3) = 0.0409 × 重量比濃度 (ppm) × 分子量
    def output_fn(gastype,mc_value):
        if gastype is None: return None
        if gastype not in MOLECULAR_MASS: return None
        return (0.0409*mc_value*input_unit_const*MOLECULAR_MASS[gastype]) *output_unit_const
    return output_fn
def vmdToMc(input_unit_const,output_unit_const):
    #公式：重量比濃度 (ppm) = 24.45 × 重量體積密度 (mg/m3) ÷ 分子量
    def output_fn(gastype,vmd_value):
        if gastype is None: return None
        if gastype not in MOLECULAR_MASS: return None
        return (24.45*vmd_value*input_unit_const / MOLECULAR_MASS[gastype]) *output_unit_const
    return output_fn

def u1ToU2(output_unit_const):
    def output_fn(gastype,vmd_value):
        return vmd_value*output_unit_const
    return output_fn

CONCENTRATION_UNIT_CONVERSION_FN={
    'vmd_ugm3':{
        'mc_ppb':vmdToMc(0.001,1000),
        'mc_ppm':vmdToMc(0.001,1),
        'vmd_mgm3':u1ToU2(0.001),
        'pm_mgm3':u1ToU2(0.001),
        'vmd_ugm3':u1ToU2(1),
        'pm_ugm3':u1ToU2(1),
    },
    'vmd_mgm3':{
        'mc_ppb':vmdToMc(1,1000),
        'mc_ppm':vmdToMc(1,1),
        'vmd_mgm3':u1ToU2(1),
        'pm_mgm3':u1ToU2(1),
        'vmd_ugm3':u1ToU2(1000),
        'pm_ugm3':u1ToU2(1000),
    },
    'mc_ppb':{
        'mc_ppb':u1ToU2(1),
        'mc_ppm':u1ToU2(0.001),
        'vmd_mgm3':mcToVmd(0.001,1),
        'pm_mgm3':mcToVmd(0.001,1),
        'vmd_ugm3':mcToVmd(0.001,1000),
        'pm_ugm3':mcToVmd(0.001,1000),
    },
    'mc_ppm':{
        'mc_ppb':u1ToU2(1000),
        'mc_ppm':u1ToU2(1),
        'vmd_mgm3':mcToVmd(1,1),
        'pm_mgm3':mcToVmd(1,1),
        'vmd_ugm3':mcToVmd(1,1000),
        'pm_ugm3':mcToVmd(1,1000),
    },
}


EC_UNIT_CONVERSION_FN={
    "ec_ampe":{
        "ec_ampe":1,
        "ec_a":1000,
        "ec_ma":1000,
        "ec_ua":1000000,
        "ec_na":1000000000,
    },
    "ec_a":{  # a = ma
        "ec_ampe":0.001,
        "ec_a":1,
        "ec_ma":1,
        "ec_ua":1000,
        "ec_na":1000000,
    },
    "ec_ma":{
        "ec_ampe":0.001,
        "ec_a":1,
        "ec_ma":1,
        "ec_ua":1000,
        "ec_na":1000000,
    },
    "ec_ua":{
        "ec_ampe":0.000001,
        "ec_a":0.001,
        "ec_ma":0.001,
        "ec_ua":1,
        "ec_na":1000,
    },
    "ec_na":{
        "ec_ampe":0.000000001,
        "ec_a":0.000001,
        "ec_ma":0.000001,
        "ec_ua":0.001,
        "ec_na":1,
    }
}

def unit_converter(input_unit,input_value,output_unit,gas_type=None):
    if input_unit not in CONCENTRATION_UNIT_CONVERSION_FN: return None
    if output_unit not in CONCENTRATION_UNIT_CONVERSION_FN[input_unit]: return None
    return CONCENTRATION_UNIT_CONVERSION_FN[input_unit][output_unit](gas_type,input_value)

def unit_converter_ppb(input_unit,input_value,gas_type=None):
    if input_unit not in CONCENTRATION_UNIT_CONVERSION_FN: return None
    return CONCENTRATION_UNIT_CONVERSION_FN[input_unit]['mc_ppb'](gas_type,input_value)

def unit_converter_na(input_unit,input_value):
    if input_unit not in EC_UNIT_CONVERSION_FN: return None
    return EC_UNIT_CONVERSION_FN[input_unit]['ec_na']*input_value




class RefSurveiesStorage:
    AWS_DYNAMODB_TABLE_NAME = 'tabairbox-dev'
    AWS_REGION_NAME = "ap-southeast-1"
    DROP_COLS=['disposed','push_time','puuid','p_time','lat','lon','setting']
    def __init__(self,key_path,ref_surveies_pool_path):
        with open(key_path, "r", encoding='UTF-8') as f:
            rows = list(csv.reader(f))[1][2:4]
            self.AWS_ACCESS_KEY_ID= rows[0]
            self.AWS_SECRET_ACCESS_KEY= rows[1]
            self.ref_surveies_pool_path=ref_surveies_pool_path
            
        class RawpushDynamoModel(dynamoModel):
            puuid = UnicodeAttribute(hash_key=True) #snUuid_y_m_ver
            m_time = NumberAttribute(range_key=True)
            p_time = NumberAttribute()
            sn= UnicodeAttribute()
            lat = NumberAttribute(default=None,null=True) 
            lon = NumberAttribute(default=None,null=True)  
            measure_time = UnicodeAttribute()
            push_time = UnicodeAttribute()
            setting = JSONAttribute(default=None,null=True)
            disposed= BooleanAttribute(default=False)
            sensors=ListAttribute()  # save as dict[] in Dynamodb,  and pull prase as json
            class Meta:
                table_name = self.AWS_DYNAMODB_TABLE_NAME
                region = self.AWS_REGION_NAME
                aws_access_key_id = self.AWS_ACCESS_KEY_ID
                aws_secret_access_key = self.AWS_SECRET_ACCESS_KEY
        self.STORE=RawpushDynamoModel
    def dyHash(self,m_time:datetime.datetime,snUuidStr:str)->str:
        return "{2}_{0}-{1:0>2d}_v0".format(m_time.year,m_time.month,snUuidStr)
    def get_HEL(self,device_uuid,m_earlier:datetime.datetime,m_later:datetime.datetime):
        dyHash=self.dyHash(m_earlier,device_uuid,)
        m_float_earlier=any_datetime_2_utc_timestamp(m_earlier)
        m_float_later=any_datetime_2_utc_timestamp(m_later)
        dir_path=os.path.join(self.ref_surveies_pool_path,dyHash)
        csv_path=os.path.join(dir_path,f"{str(m_float_earlier)}.csv")
        return {
            'hash':dyHash,
            'm_later': m_float_later,
            'm_earlier':m_float_earlier,
            'dir_path':dir_path,
            'csv_path':csv_path
        }
    def downloadDataOneDay(self,hel)->pd.DataFrame:
        # print(self.dyHash(m_earlier,device_uuid,),)
        # print(any_datetime_2_utc_timestamp(m_earlier),any_datetime_2_utc_timestamp(m_later))
        # return
        m=utc_timestamp_2_utc_datetime(hel['m_earlier'])
        print(f'download {m} {hel["csv_path"]}')
        q= self.STORE.query(
                hel['hash'],
                range_key_condition=(
                    self.STORE.m_time.between(hel['m_earlier'],hel['m_later'])
                    ),
                limit=5000
                )
        def to_json(d):
            d["setting"]=json.dumps(d["setting"],ensure_ascii=False)
            d["sensors"]=json.dumps(d["sensors"],ensure_ascii=False)
            return d
        df=pd.DataFrame([ to_json(item.attribute_values) for item in q ])
        if not os.path.exists(hel['dir_path']):
            os.makedirs(hel['dir_path'])
        df.to_csv(hel['csv_path'],index=False,encoding='utf-8')
        return df
        
    def getData(self,device_uuid,now:datetime.datetime,m_earlier:datetime.datetime,m_later:datetime.datetime)->pd.DataFrame:
        m_start_chunk=m_earlier.replace(hour=0,minute=0,second=0,microsecond=0)
        iday=-1
        latest_chunk=now.replace(hour=0,minute=0,second=0,microsecond=0)
        datas=[]
        while True:
            iday+=1
            this_chunk=m_start_chunk+relativedelta(days=iday)
            end_chunk=this_chunk+relativedelta(days=1)+relativedelta(seconds=-1)
            if this_chunk >=m_later:
                break
            hel=self.get_HEL(device_uuid,this_chunk,end_chunk)
            if self.checkIsStoredLocal(hel) and latest_chunk != this_chunk :
                datas.append(pd.read_csv(hel['csv_path']))
            else:
                datas.append(self.downloadDataOneDay(hel))
        df=pd.concat(datas,axis=0, ignore_index=True)
        df=self.df_process( df)
        mask = (df['m_time'] > any_datetime_2_utc_timestamp(m_earlier) ) & (df['m_time'] < any_datetime_2_utc_timestamp(m_later) )
        df=df[mask]
        df.reset_index(drop=True, inplace=True)
        return df
    def df_process(self, df):
        #print(df)
        def serializer(x):
            l=json.loads(x)
            d={}
            for item in l:
                d[item["name"]]=item["value"]
            return d
        df2=pd.DataFrame(df["sensors"].map(serializer).tolist())
        df.drop(self.DROP_COLS+["sensors"], axis=1, inplace=True)
        return pd.concat([df,df2], join='outer',axis=1)
    def checkIsStoredLocal(self,hel):
        return os.path.exists(hel['csv_path'])



  



class RefStationsStorage:
    latest_name='resource'
    reg_yyyy_mm='(19|20)\d\d[-](0[1-9]|1[012])'
    reg_yyyy_mm_dd='(19|20)\d\d[-](0[1-9]|1[012])[-](0[1-9]|[12][0-9]|3[01])'
    tz=pytz.timezone('Asia/Taipei')
    def __init__(self,ref_stations_pool_path):
        self.ref_stations_pool_path=ref_stations_pool_path
        
    def isDirInZip(self,zfio,file_path):
        file_path=file_path.replace('\\','/')
        return any(x.startswith("%s/" % file_path.rstrip("/")) for x in zfio.namelist())
    def isFileInZip(self,zfio,file_path):
        file_path=file_path.replace('\\','/')
        return file_path in zfio.namelist()
    def zipIOStrLoad(self,zfio,file_path):
        file_path=file_path.replace('\\','/')
        with zfio.open(file_path) as f:
            return f.read()

    def zipIOLoad(self,zfio,file_path):
        file_path=file_path.replace('\\','/')
        return zfio.open(file_path)

    def df_process(self,df,sort=False):
        def serializer(x):
            d={}
            for i,name in enumerate( x['itemengname']):
                try:
                    v=unit_converter_ppb("mc_"+x["itemunit"][i],float(x["concentration"][i]),gas_type=name)
                except:
                    v=np.NaN
                d[ "{}#{}".format(name.lower(),'ppb')]=v 
            return d
        df.columns = map(str.lower, df.columns)
        df=df.groupby('monitordate').agg({'itemid':list,'itemname':list,'itemengname':list,'itemunit':list,'concentration':list,'siteid':'first','sitename':'first','county':'first'}).reset_index()
        dcolumns = ['itemid', 'itemname','itemengname', 'itemunit', 'concentration']
        df2=pd.DataFrame(pd.Series( df[dcolumns].to_dict(orient='records')).map(serializer).tolist())
        df.drop(dcolumns, axis=1, inplace=True)
        df=pd.concat([df,df2], join='outer',axis=1)
        df['monitordate']=pd.to_datetime(df['monitordate'], format='%Y-%m-%d %H:%M:%S',infer_datetime_format=True).dt.tz_localize(self.tz).dt.tz_convert('UTC')
        # if sort: # no need
        #     df.sort_values('monitordate',axis=0,inplace=True)
        #print(df)
        return df

    def get_zip_target_from_latest_file(self,station_name,chunk_date_str,this_chunk,latest_chunk)->pd.DataFrame:
        # first search AQX_P_190_Resource 
        chunk_file="{}_{}.zip".format(station_name,self.latest_name)
        zip_path=os.path.join(self.ref_stations_pool_path,chunk_file)
        if os.path.isfile(zip_path):
            zfio=zipfile.ZipFile(zip_path)
            print(latest_chunk, this_chunk)
            if latest_chunk == this_chunk:
                for scan_target in zfio.namelist():
                    if not re.search(self.reg_yyyy_mm, scan_target):
                        print( f"load zip: {zip_path} - {scan_target}")
                        return self.df_process(pd.read_csv(self.zipIOLoad(zfio,scan_target)))
                    
            else:
                for scan_target in zfio.namelist():
                    if chunk_date_str in scan_target:
                        print( f"load zip: {zip_path} - {scan_target}")
                        return self.df_process(pd.read_csv(self.zipIOLoad(zfio,scan_target)))
        return None

    def get_zip_target_from_date_file(self,station_name,chunk_date_str,this_chunk,latest_chunk)->pd.DataFrame:
        # 2nd search AQX_P_190_20xx_xx
        chunk_file="{}_{}.zip".format(station_name,chunk_date_str)
        zip_path=os.path.join(self.ref_stations_pool_path,chunk_file)
        if os.path.isfile(zip_path):
            zfio=zipfile.ZipFile(zip_path)
            namelist=sorted(zfio.namelist())
            for scan_target in namelist:
                if not re.search(self.reg_yyyy_mm_dd, scan_target):
                    print( f"load zip: {zip_path} - {scan_target}")
                    return self.df_process(pd.read_csv(self.zipIOLoad(zfio,scan_target)))
            pds=[]
            for scan_target in namelist:
                searched=re.search(self.reg_yyyy_mm_dd, scan_target)
                if searched:
                    #print(searched.group(0),this_chunk)
                    print( f"load zip: {zip_path} - {scan_target}")
                    pds.append( self.df_process(pd.read_csv(self.zipIOLoad(zfio,scan_target)),sort=True))
            return pd.concat(pds,axis=0, ignore_index=True)
        return None

    def getData(self,station_name,now:datetime.datetime,m_earlier:datetime.datetime,m_later:datetime.datetime)->pd.DataFrame:
        m_start_chunk=m_earlier.astimezone(self.tz).replace(day=1,hour=0,minute=0,second=0,microsecond=0)
        latest_chunk=now.astimezone(self.tz).replace(day=1,hour=0,minute=0,second=0,microsecond=0)
        imonth=-1
        datas=[]
        
        while True:
            imonth+=1
            this_chunk=m_start_chunk+relativedelta(months=imonth)
            chunk_date_str="{0}-{1:0>2d}".format(this_chunk.year,this_chunk.month)
            #end_chunk=this_chunk+relativedelta(months=1)+relativedelta(hours=-1)
            if this_chunk >=m_later:
                break
            zip_df=self.get_zip_target_from_latest_file(station_name,chunk_date_str,this_chunk,latest_chunk)
            if zip_df is None:
                zip_df=self.get_zip_target_from_date_file(station_name,chunk_date_str,this_chunk,latest_chunk)
            if zip_df is None:
                raise Exception(f"{self.ref_stations_pool_path}, you miss the date of {chunk_date_str}")
            datas.append(zip_df)
        df=pd.concat(datas,axis=0, ignore_index=True)
        mask = (df['monitordate'] >= m_earlier ) & (df['monitordate'] <= m_later )
        df=df[mask]
        df.reset_index(drop=True, inplace=True)
        return df
        

class RefSurveiesVender:
    def __init__(self,surv_storage:RefSurveiesStorage,device_uuid,now,m_earlier,m_later):
        self.datas=surv_storage.getData(device_uuid,now,m_earlier,m_later)
        self.now=now
        self.device_uuid=device_uuid
        self.m_earlier=m_earlier
        self.m_later=m_later
        
    def getData(self):
        return self.datas

    def getDataByDateList(self,regression_date_list)->List [pd.DataFrame]:
        """shift back 1~2 hour to o'clock  with whole 1 hour data
        """

        r=[]
        for date in regression_date_list:
            new_date_end=date.replace(minute=0,second=0,microsecond=0)+relativedelta(hours=-1)
            new_date_start=date.replace(minute=0,second=0,microsecond=0)+relativedelta(hours=-2)
            if self.m_earlier > new_date_start:
                raise Exception(f"request the earliest date is {new_date_start}")
            mask = (self.datas['m_time'] >= any_datetime_2_utc_timestamp(new_date_start) ) & (self.datas['m_time'] <= any_datetime_2_utc_timestamp(new_date_end) )
            r.append(self.datas[ mask  ])
        return r
    def get_value_array_of_na(self,sensor_order,df:pd.DataFrame)->np.ndarray:
        def get(x):
            return x.values[0]
        if isinstance(df,pd.Series):
            def get(x):
                return x
        def fn(name):
            unittype=get_sesnor_unit(name)
            input=get(df[name])
            if pd.isna(input) :
                raise MissingValueException
            return unit_converter_na(unittype,input) 
        return np.asarray([ fn(i)  for i in sensor_order ])
    
    def check_diff_zreo(self,point_surv_i:np.ndarray,point_surv_j:np.ndarray)->bool:
        return np.any(np.isin(point_surv_i-point_surv_j,[0]))

class RefStationsVender:
    def __init__(self,stat_storage:RefSurveiesStorage,station_name,now,m_earlier,m_later):
        self.datas=stat_storage.getData(station_name,now,m_earlier,m_later)
        self.now=now
        self.station_name=station_name
        self.m_earlier=m_earlier
        self.m_later=m_later
    def getData(self):
        return self.datas

    def getDataByDateList(self,regression_date_list)->List [pd.DataFrame]:
        """shift back 1~2 hour to o'clock
        """
        r=[]
        for date in regression_date_list:
            new_date=date.replace(minute=0,second=0,microsecond=0)+relativedelta(hours=-1)
            if self.m_earlier > new_date:
                raise Exception(f"request the earliest date is {new_date}")
            r.append(self.datas.loc[self.datas['monitordate'] == new_date])
        return r
    def get_value_array_of_ppb(self,gas_order,df:pd.DataFrame)->np.ndarray:
        def get(x):
            return x.values[0]
        if isinstance(df,pd.Series):
            def get(x):
                return x
        def fn(name):
            for h in df.columns:
                col_name_s=h.split('#')
                if len(col_name_s) <=1 : continue
                if col_name_s[0] == name:
                    unittype= 'mc_'+col_name_s[1]
                    input=get(df[h])
                    if pd.isna(input):
                        raise MissingValueException
                    return unit_converter_ppb(unittype,input,gas_type=name)   
        return np.asarray([ fn(i)  for i in gas_order ])
       


def getIndex_of_current_date_and_before(current_date:Union[datetime.datetime,str],date_list):
    date=None
    if isinstance(current_date,datetime.datetime):
        date=current_date
    elif isinstance(current_date,str):
        date=any_isoformat_2_utc_datetime(current_date) 
    else:
        raise Exception(f"wrong current_date type: {type(current_date)}")

    l=[ i for i,item in enumerate( date_list) if item <= current_date]
    return l
            
        

