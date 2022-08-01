import os 
import datetime
import pandas as pd
from typing import List
import matplotlib.pyplot as plt  
import matplotlib.dates as mdate
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

GAS_ORDER=["co","so2","o3","no2"]
SENSOR_ORDER=[
            "c1_0#ec_na#0",
            "c1_5#ec_na#0",
            "c1_4#ec_na#0",
            "c1_3#ec_na#0"
            ]

COLOR_MAP = ['#FFBE7D', '#4E79A7', '#F28E2B', '#A0CBE8', '#59A14F', '#8CD17D', '#B6992D', '#F1CE63', '#499894', '#86BCB6', '#E15759', '#FF9D9A', '#79706E', '#BAB0AC', '#D37295', '#FABFD2', '#B07AA1', '#D4A6C8', '#9D7660', '#D7B5A6']
FIG_PATH = ['fig', 'fig\\station', 'fig\\device_current', 'fig\\device_concentration']

SENSOR_CODE_MAP = {
    'co':{'station':'co#ppb', 'device_concentration':'pred_co#ppb', 'device_current':'c1_0#ec_na#0'},
    'so2':{'station':'so2#ppb', 'device_concentration':'pred_so2#ppb', 'device_current':'c1_5#ec_na#0'},
    'o3':{'station':'o3#ppb', 'device_concentration':'pred_o3#ppb', 'device_current':'c1_4#ec_na#0'},
    'no2':{'station':'no2#ppb', 'device_concentration':'pred_no2#ppb', 'device_current':'c1_3#ec_na#0'}
}


def make_fig_folder():
    pass


class PlotStation:
    make_fig_folder()
    def __init__(self,stat_df:pd.DataFrame):
        self.df = stat_df
        self.time_cols = 'monitordate'
        self.fig_path = 'fig\\station'
        self.target = self.fig_path.split('\\')[-1]
        self.target_cols=GAS_ORDER
        self.interval = 5
        self.unit = 'ppb'
        self.ylim_flag = False
        pass
    
    def process_date(self,date_cols):
        return date_cols
    
    def plot_single(self):
        '''
        individual fig for each sensor
        '''
        for index, value in enumerate(self.target_cols):
            fig,ax = plt.subplots(figsize = (8,6))
            plt.plot(self.process_date(self.df[self.time_cols]), self.df[SENSOR_CODE_MAP[value][self.target]], label = SENSOR_CODE_MAP[value][self.target], color = COLOR_MAP[index])
            plt.xlabel(self.time_cols, fontsize=14, fontweight='bold')
            plt.ylabel(SENSOR_CODE_MAP[value][self.target], fontsize=14, fontweight='bold')
            plt.legend(bbox_to_anchor=(1.2,0.1), loc="lower right")
            plt.grid()
            if value == 'co' and self.ylim_flag:
                plt.ylim(0, 1000)
            plt.savefig(f"{self.fig_path}\\Single_{value}_{self.unit}.png", bbox_inches='tight')
            plt.clf()
    
    
    
    
    
    
    
    
    
    
class PlotDeviceConcentration(PlotStation):
    def __init__(self, full_df):
        super().__init__(full_df)
        self.df = full_df
        self.time_cols = 'measure_time'
        self.fig_path = 'fig\\device_concentration'
        self.target = self.fig_path.split('\\')[-1]
        self.interval = 150
        self.unit = 'ppb'
        self.ylim_flag = True
        self.temp_colsname = 'c0_tm0#tm_c#0'
        self.rh_colsname = 'c0_rh0#rh_p#0'
        
    def process_date(self,date_cols):
        return pd.to_datetime(date_cols)
    
    # def plot_temp(self):
    #     print(self.df.columns)

    #     fig,ax = plt.subplots(figsize = (8,6))
    #     plt.plot(self.process_date(self.df[self.time_cols]), self.df[self.temp_colsname], label = "Temp", color = 'orange')
    #     plt.xlabel(self.time_cols, fontsize=14, fontweight='bold')
    #     plt.legend(bbox_to_anchor=(1.2,0.1), loc="lower right")
    #     plt.grid()
    #     plt.show()
    #     # plt.savefig(f"{self.fig_path}\\Single_{value}_{self.unit}.png", bbox_inches='tight')
    #     plt.clf()
    
    

            
class PlotDeviceCurrent(PlotDeviceConcentration):
    def __init__(self, full_df):
        super().__init__(full_df)
        self.fig_path = 'fig\\device_current'
        self.target = self.fig_path.split('\\')[-1]
        self.unit = 'nA'
        self.ylim_flag = False
