import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr
import math
from scipy.stats import gaussian_kde
from matplotlib import cm
import operator
import csv

#   ____ _                  ___     _   _      _                     
#  / ___| | __ _ ___ ___   ( _ )   | | | | ___| |_ __   ___ _ __ ___ 
# | |   | |/ _` / __/ __|  / _ \/\ | |_| |/ _ \ | '_ \ / _ \ '__/ __|
# | |___| | (_| \__ \__ \ | (_>  < |  _  |  __/ | |_) |  __/ |  \__ \
#  \____|_|\__,_|___/___/  \___/\/ |_| |_|\___|_| .__/ \___|_|  |___/
#                                               |_|                  

class Gene:

    def __init__(self, data):
    	self.name = data[0]
    	valueList = data[2:]
    	valueList = [float(i) for i in valueList ]
    	self.dis_mean = (valueList[0]+valueList[1]+valueList[2])/3
    	self.point_five_mean = (valueList[3]+valueList[4]+valueList[5])/3
    	self.two_mean = (valueList[6]+valueList[7]+valueList[8])/3
    	self.mat_mean = (valueList[9]+valueList[10]+valueList[11])/3
    	
    def getMeanList(self):
        return [self.dis_mean, self.point_five_mean, self.two_mean, self.mat_mean]
        
    def calCorrelation(self, refSample):
        corr = pearsonr(self.getMeanList(),refSample)
        self.r_val = corr[0]
        self.p_val = corr[1]
        
        
    def get_r_val(self):
        return self.r_val
        
    def get_p_val(self):	
        return self.p_val
        
    def isPositiveTrend(self, threshold):
        if(math.isnan(self.r_val)):
            return False
        
        if(self.r_val > threshold and self.p_val < 0.05):
            return True
            
        return False
        
    def isNegativeTrend(self, threshold):
        if(math.isnan(self.r_val)):
            return False
        
        if(self.r_val < threshold and self.p_val < 0.05):
            return True
            
        return False
        
    def __lt__(self, other):
        return self.r_val < other.r_val


def onePlot(yAxis, color = "red", is_show = True, name="spl9", is_save=False):
    x = [ "disected", "0.5", "2", "mature"]
    xnum = [1, 2, 3, 4]
    plt.scatter(x, yAxis, color=color)
    z = np.polyfit(xnum, yAxis, 1)
    p = np.poly1d(z)
    plt.plot(x, yAxis, "r--", label="spl9")
    if(is_show):
        plt.show()
        
    if(is_save):
        name = "plots/" + name + ".png"
        plt.savefig(name, bbox_inches='tight')
    plt.cla()
    

def getCorrelation(geneSample, refSample):
    return pearsonr(geneSample,refSample)


def kde(density_list):
    density = gaussian_kde(density_list) # x: list of price
    xgrid = np.linspace(min(density_list), max(density_list), len(density_list))   
    plt.plot(xgrid, density(xgrid))
    plt.show()
    
def multiPlot(trend_list, title = "Positive Trend", is_save=False, ref_gene = None):
        
    if(title == "Positive Trend"):
        trend_list = sorted(trend_list, key=operator.attrgetter('r_val'), reverse=True)
    else:
        trend_list = sorted(trend_list, key=operator.attrgetter('r_val'))

    x = [ "disected", "0.5", "2", "mature"]
    xnum = [1, 2, 3, 4]
    
    length = len(trend_list)
    if(length>10):
        length = 10
        
    if(ref_gene):
        trend_list.append(ref_gene)
        length = length+1

    for index in range(length):
        col = (np.random.random(), np.random.random(), np.random.random())
        plt.scatter(x, trend_list[index].getMeanList(), col)
        z = np.polyfit(xnum, trend_list[index].getMeanList(), 1)
        p = np.poly1d(z)
        label = trend_list[index].name
        if(label == reference_data):
            label = "SPL9"
        plt.plot(x, trend_list[index].getMeanList(), label=label)
        
    plt.title(title)
    plt.legend(loc="upper right")
    plt.xlabel("Types")
    plt.ylabel("Tpms")
    
    if(is_save):
        name = "data/" + title + ".png"
        plt.savefig(name, bbox_inches='tight')
    plt.show()
    plt.cla()




#   ____                       _        _   _     _   _           _                                        _     _                               _         _       _   _   _             
#  / ___| ___ _ __   ___   ___| |_ __ _| |_(_)___| |_(_) ___ __ _| |   ___ ___  _ __ ___  _ __   __ _ _ __(_)___(_) ___  _ __     __ _ _ __   __| |  _ __ | | ___ | |_| |_(_)_ __   __ _ 
# | |  _ / _ \ '_ \ / _ \ / __| __/ _` | __| / __| __| |/ __/ _` | |  / __/ _ \| '_ ` _ \| '_ \ / _` | '__| / __| |/ _ \| '_ \   / _` | '_ \ / _` | | '_ \| |/ _ \| __| __| | '_ \ / _` |
# | |_| |  __/ | | |  __/ \__ \ || (_| | |_| \__ \ |_| | (_| (_| | | | (_| (_) | | | | | | |_) | (_| | |  | \__ \ | (_) | | | | | (_| | | | | (_| | | |_) | | (_) | |_| |_| | | | | (_| |
#  \____|\___|_| |_|\___| |___/\__\__,_|\__|_|___/\__|_|\___\__,_|_|  \___\___/|_| |_| |_| .__/ \__,_|_|  |_|___/_|\___/|_| |_|  \__,_|_| |_|\__,_| | .__/|_|\___/ \__|\__|_|_| |_|\__, |
#                                                                                        |_|                                                        |_|                            |___/ 



min_threshold = -0.75
max_threshold = 0.75
reference_data = "CARHR137240" # spl9 gene

df = pd.read_csv("data.csv")
valueList = df.values.tolist()
valueListNumpy = np.array(valueList)
dataList = valueListNumpy
dataMap = {}
for row in dataList:
    dataMap[row[0]] = Gene(row)

# plotting spl9 trend
refList = dataMap[reference_data].getMeanList();
#onePlot(refList)


# Find correlation for all
r_val_list = []
p_val_list = []
for key, val in dataMap.items():
    val.calCorrelation(refList)
    r_val_list.append(val.get_r_val())
    p_val_list.append(val.get_p_val())
    

#for_kde = [i for i in r_val_list if not math.isnan(i)]
#kde(for_kde)

positive_trend_list = []
negative_trend_list = []

for key, val in dataMap.items():
    if(val.isPositiveTrend(max_threshold)):
        positive_trend_list.append(val)
    elif( val.isNegativeTrend(min_threshold)):
        negative_trend_list.append(val)
 
print("=======================================================================================")       
data = [i for i in positive_trend_list]
print("Positive trend numbers : {}".format(len(data)))

data_Negative = [i for i in negative_trend_list]
print("Negative trend numbers : {}".format(len(data_Negative)))
print("=======================================================================================")       

multiPlot(positive_trend_list, title = "Positive Trend", is_save=True, ref_gene=dataMap[reference_data])
multiPlot(negative_trend_list, title = "Negative Trend", is_save=True, ref_gene=dataMap[reference_data])


# Save to csv file

#Positive trend
fileName = "data/trend_" + str(max_threshold) + ".csv"
with open(fileName, 'w', encoding='UTF8') as f:
    writer = csv.writer(f)
    writer.writerow(["GeneID", "R-Value", "P-Value", "Correlation-Type"])
    for gene in positive_trend_list:
        writer.writerow([gene.name, gene.get_r_val(), gene.get_p_val(), "Positive"])
    for gene in negative_trend_list:
        writer.writerow([gene.name, gene.get_r_val(), gene.get_p_val(), "Negative"])
