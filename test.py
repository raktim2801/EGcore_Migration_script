#!/usr/bin/env python
# coding: utf-8

"""
A script to automatize Studies from `EGcore` 
into template suitable for Onboarding in Analysis hub.
"""

# ## Load modules
import pandas as pd
#import csv
import numpy as np
import logging
import json
#import warnings
#from pathlib import Path
import os
import subprocess
import sys
import argparse
import time
import io
import shutil


# ## Get Study ID
parser = argparse.ArgumentParser(
    description= "A script to automatize Studies from `EGcore` into template suitable for Onboarding in Analysis hub."
)
parser.add_argument('-i', '--id', required=True, help='The source Eaglecore Study ID')
parser.add_argument('-e', '--errors', required=False, help='The mode of handleing wrong DataType values, default= raise', 
                    choices=['coerce', 'raise', 'ignore'], default='raise')
parser.add_argument('-d', '--debug', help="Print lots of debugging statements", 
                    action="store_const", dest="loglevel", const=logging.DEBUG, default=logging.WARNING)
parser.add_argument('-v', '--verbose', help="Be verbose",
                    action="store_const", dest="loglevel", const=logging.INFO)
parser.add_argument('-dl', '--delim', required=False, help='The delimiter to be used, default= False',
                    default=False)
args = parser.parse_args()

SID = args.id
errors = args.errors
delim=args.delim

print("----------------------------------------------------------------------------------------------")
print(SID)
print("Level of logging: ", args.loglevel)


timestr = time.strftime("%Y-%m-%d_%H_%M_%S")
dir=f'./validationReports/egcore{SID}/log'
try:
    os.makedirs(dir)
    print(f"Folder {dir} created.")
except Exception:
    print("Folder already exists.")
logging.basicConfig(filename=f'./validationReports/egcore{SID}/log/{SID}_{timestr}.log', 
                    format='%(levelname)s>> %(asctime)s: %(message)s', 
                    datefmt='%d/%m/%Y %I:%M:%S %p',
                    encoding='utf-8', level=args.loglevel)

logging.info(f'STARTED at: {time.strftime("%a, %d %b %Y %H:%M:%S")}')


buf = io.StringIO()


############################################################################################

def CREATE_VALIDATION_DICTS():
    """# ## Load Dictionaries input file

    return: ddsyn2, ddSIM, ddEntity, ddDType
    """
    try:
        fp="/Users/raktimmaiti_mbp/Documents/GitHub/EagleGenomics-CV/master/eaglegenomics-cv.json"
        with open(fp) as f:
            cv=json.load(f)
        logging.info('Input mapping file successfully loaded')
    except:
        logging.error("""CV file not loaded""")
        raise ImportError('CV file not loaded')

    """Create Dicts from input file"""
    """# ### Load Synonym dict """
    ddsyn2={}
    for x in cv['classes'][1]['terms'].keys():
        temp=cv['classes'][1]['terms'][x]["synonyms"].strip().split(",")
        temp=[y.strip() for y in temp if len(y) > 0]
        temp=[y.strip() for y in temp if len(y) > 0]
        for zz in temp:
            ddsyn2[zz] = x

    """# ### Load Entity and SIM mapping"""
    ddEntity={x:cv['classes'][1]['terms'][x]['entity'] for x in cv['classes'][1]['terms'].keys()}
    ddSIM={x:cv['classes'][1]['terms'][x]['default_display'] for x in cv['classes'][1]['terms'].keys()}
    ddDType={x:cv['classes'][1]['terms'][x]['datatype'] for x in cv['classes'][1]['terms'].keys()}
    
    """# ### Load valid values mapping"""
    ddNumRange={}
    ddCatRange={}
    for x in cv['classes'][1]['terms'].keys():
        y=cv['classes'][1]['terms'][x]['valid_values']
        z=cv['classes'][1]['terms'][x]['datatype']
        if not pd.isna(y):
            #print(x,cv['classes'][1]['terms'][x]['valid_values'])
            if z == 'Numeric':
                ddNumRange[x]=y
            else:
                ddCatRange[x]=y
                
    ddNumRange={x:eval(ddNumRange[x]) for x in ddNumRange}
    ddCatRange={x:[_.strip() for _ in ddCatRange[x].split(',')] for x in ddCatRange}
        
    return ddsyn2, ddSIM, ddEntity, ddDType, ddNumRange, ddCatRange


"""# ## Download data from S3 bucket"""
def createMappingURL(sid, outpath="./"):
    """Create mapping.csv file urls from study id and the output file path"""
    url = f"aws s3 cp s3://eaglecore-unilever/studies/{sid}/mapping.csv {outpath}{sid}_mapping.csv"
    return url

def createStudyDetailsURL(sid, outpath="./"):
    """Create StudyDetails.json file urls from study id and the output file path"""
    url = f"aws s3 cp s3://eaglecore-unilever/studies/{sid}/StudyDetails.json {outpath}{sid}_StudyDetails.json"
    return url

"""# ### Load mapping.csv file func"""
def LOAD_MAPPING_FILE(SID, ddsyn2):
    """# ### Load mapping.csv file into a pandas DataFrame.
    process the DataFrame and map the column headers.
    Remove duplicate columns
    Drops Blank columns.

    return: a DataFrame
    """
    x=f"./validationReports/egcore{SID}/{SID}_mapping.csv"
    delimit=[r'[,;\t]', ",", "\t"]
    logging.info("Starting to load `mapping.csv` file.")
    
    try:
        with open(x) as f:
            ID= SID
            
            if delim:
                dd=pd.read_csv(f, delimiter=delim, engine='python', dtype=object)
            else:
                for fn in delimit:
                    try:
                        dd=pd.read_csv(f, delimiter=fn, engine='python', dtype=object)
                        break
                    except TypeError:
                        logging.warning(f"{ID} > csv not loaded, tried {fn}")
                        continue
            #dd=pd.read_csv(f, delimiter=r'[,;\t]', engine='python')
            logging.info("`mapping.csv` file SUCCESSFULLY loaded.")

            cc=list(dd.columns)
            cc2=[]
            for ccx in cc:
                if ccx in ddsyn2.keys():
                    dx=ddsyn2[ccx]
                    while dx in cc2:
                        dx=dx+'_|'
                    cc2.append(dx)
                else:
                    cc2.append("~"+ ccx)
            dd.columns=cc2
    except FileNotFoundError:
        logging.critical("NO mapping file found in EGCORE.")
        sys.exit("NO mapping file found in EGCORE.")

    #dd.insert(loc=0, column='#EGcoreStudyID', value=ID)
    dd.reset_index(inplace=True, drop=True)
    cd=dd.loc[:,dd.columns.duplicated()].copy()
    dd = dd.loc[:,~dd.columns.duplicated()].copy()
    mappingFile= dd.copy()
    logging.info("DataFrame created. Attributes Mapped and duplicates removed.")
    logging.info(f"DataFrame size: {mappingFile.shape}")
    print(f"DataFrame size: {mappingFile.shape}")

    # ### Delete Blank columns
    mappingFile = mappingFile.dropna(axis=1, how='all')
    logging.info("DataFrame: BLANK Attributes removed.")
    logging.info(f"DataFrame new size: {mappingFile.shape}")
    
    # ### strip all values
    mappingFile = mappingFile.apply(lambda x: x.apply(lambda y: y.strip() if type(y) == type('') else y), axis=0)
    return mappingFile


def convertDtype(dseries, ddDType, errors=errors):
    s=dseries.name
    if s in ddDType:
        dt = ddDType[s]
    else:
        dt = "String"
    
    if dseries.dtype == 'float64' or dseries.dtype == 'int64':
        dt2 = "Numeric"
    else:
        dt2 = "String"
    
    #if dt != dt2:
        #logging.warning(f"For '{s}' column, expected dataType {dt}, actual dataType {dt2}.")
    
    if dt == 'DateTime':
        dseries = dseries.apply(lambda x: pd.to_datetime(x, format='%d_%m_%Y', errors=errors))
        logging.info(f"converted dataType of {s}: {dt2} --> DateTime, error mode = {errors}")
        #dseries = dseries.apply(lambda x: x.strftime('%d-%b-%Y'))
    elif dt == 'Numeric':
        dseries = dseries.apply(lambda x: pd.to_numeric(x, errors=errors))
        logging.info(f"converted dataType of {s}: {dt2} --> Numeric, error mode = {errors}")
    else:
        #dseries = dseries.fillna('').apply(str)
        pass
    
    return dseries


def checkValidValues(dseries, ddDType, ddNumRange, ddCatRange):
    s=dseries.name
    outOfRngCol=[]
    if s in ddDType:
        dt=ddDType[s]
    else:
        dt="String"
    
    if dt == 'Numeric' and s in ddNumRange:
        Ndata=dseries.copy()
        if Ndata.dropna().between(ddNumRange[s][0],ddNumRange[s][1]).all():
            pass
        else:
            #print(f"{s} have out of range numerical values")
            logging.error(f"Atribute {Ndata.name}: contains Out of Range Numerical Values")
            print(ddNumRange[s])
            print(Ndata)
            outOfRngCol.append(Ndata.name)
    if dt == 'String' and s in ddCatRange:
        Cdata=dseries.copy()
        if Cdata.dropna().isin(ddCatRange[s]).all():
            pass
        else:
            #print(f"{s} have out of range categorical values")
            logging.error(f"Atribute {Cdata.name}: contains Out of Range Categorical Values")
            print(ddCatRange[s])
            print(Cdata)
            outOfRngCol.append(Cdata.name)
    
    return outOfRngCol

def checkValidNumValues(dseries, ddNumRange):
    s=dseries.name
    outOfRngCol=[]
    
    if s in ddNumRange:
        Ndata=dseries.copy()
        if Ndata.dropna().between(ddNumRange[s][0],ddNumRange[s][1]).all():
            pass
        else:
            #print(f"{s} have out of range numerical values")
            logging.error(f"Atribute {Ndata.name}: contains Out of Range Numerical Values")
            print(Ndata.name, "Valid Range: ", ddNumRange[s])
            print(Ndata.name, "Actual data-> min: ", Ndata.min(),", max: ", Ndata.max())
            outOfRngCol.append(Ndata.name)
    else:
        pass
    
    return outOfRngCol


def checkValidCatValues(dseries, ddCatRange):
    s=dseries.name
    outOfRngCol=[]
    
    if s in ddCatRange:
        Cdata=dseries.copy()
        if Cdata.dropna().isin(ddCatRange[s]).all():
            pass
        else:
            #print(f"{s} have out of range categorical values")
            logging.error(f"Atribute {Cdata.name}: contains Out of Range Categorical Values")
            print(Cdata.name, "Valid Values: ", ddCatRange[s])
            #print(Cdata.name, "Actual Values: ", Cdata.unique())
            print(Cdata.name, "Unknown Values: ", set(Cdata.unique())-set(ddCatRange[s]))
            outOfRngCol.append(Cdata.name)
    else:
        pass
    
    return outOfRngCol

"""# ### Load StudyDetails.json file func"""
def LOAD_JSON_FILE(SID, ddsyn2):
    """# ### Load StudyDetails.json file.
    Creates a dictionary.

    return: the dict
    """
    x=f"./validationReports/egcore{SID}/{SID}_StudyDetails.json"
    with open(x) as infile:
        jsonFile=json.load(infile)

    temp=jsonFile.copy()
    for x in temp:
        if x in ddsyn2:
            jsonFile[ddsyn2[x]]=jsonFile.pop(x)

    return jsonFile


# ## Load BLANK Templates
def LOAD_BLANK_TP():
    """# ## Load BLANK Templates.
    returns a list of DFs
    """
    meta=pd.read_excel('./EGCORE_TP/EGCORE_CuratedMetaData.xlsx', sheet_name=None, header=None)
    metacopy=meta.copy()
    logging.info("BLANK MetaData TP loaded.")
    ### NOT LOADING Measurement TP, as currently NOT required
    logging.debug("NOT LOADING Measurement TP, as currently NOT required")
    #measure= pd.read_excel('./EGCORE_TP/curMeasureX.xlsx', sheet_name=None, header=None)
    #measurecopy=measure.copy()

    ########## DEFINE ALL SHEETS as DFs ##########
    UserGuide= meta['UserGuide']
    Config= meta['Config']
    Attributes= meta['Attributes']
    Metadata= meta['Metadata']
    mapping= meta['MAPPING']
    Relation= meta['Relation']
    study= meta['Study']
    experiment= meta['Experiment']
    treatment=meta['Treatment']
    timepoint=meta['Timepoint']
    location=meta['Location']
    subject=meta['Subject']
    sample=meta['Sample']
    measMetaX=meta['CuratedMeasurementMetaData1']
    measMetaXtra=meta['CuratedMeasurementMetaData2']

    return study, UserGuide, Config, Attributes, Metadata, Relation

print("\n--------------------------------\n")

with open("egcore_downloadfiles.sh","w") as f1:
    m=createMappingURL(SID, outpath=f"./validationReports/egcore{SID}/")
    s=createStudyDetailsURL(SID, outpath=f"./validationReports/egcore{SID}/")
    f1.write(m+'\n')
    f1.write(s)

r=subprocess.call(['sh', './egcore_downloadfiles.sh'])
print(r)
if r==0:
    logging.info("Files from S3 bucket downloaded successfully.")
else:
    logging.critical("Files from S3 bucket could not be downloaded.")
## Delete the temp file
os.remove("./egcore_downloadfiles.sh")


def main():
    dataValidation=[]
    """# ## Load Dictionaries input file"""
    ddsyn2, ddSIM, ddEntity, ddDType, ddNumRange, ddCatRange = CREATE_VALIDATION_DICTS()
    """Load input study files"""
    mappingFile = LOAD_MAPPING_FILE(SID, ddsyn2)
    jsonFile = LOAD_JSON_FILE(SID, ddsyn2)
    study, UserGuide, Config, Attributes, Metadata, Relation = LOAD_BLANK_TP()
    #study= UPDATE_STUDY_FROM_JSON(SID, study, jsonFile)
    
    # ### Check all attributes are defined or not
    allDefined=[x for x in mappingFile.columns if "~" in x or "_|" in x]
    if len(allDefined) > 0:
        logging.critical("mapping file have Undefined or Duplicate attributes.")
        print(f"Undefined or Duplicate attributes present: {allDefined}")
        print(mappingFile.head(5))
        #sys.exit("mapping.csv file have Undefined or Duplicate attributes.")
        dataValidation.append(False)
    else:
        dataValidation.append(True)
    
    # ## Mandatory attributes check (NOT IMPLEMENTED)
    dataValidation.append(True)
        
    mappingFile.info(buf=buf)
    logging.info('mapping.csv info - {}'.format(buf.getvalue()))
