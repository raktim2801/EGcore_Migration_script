import pandas as pd
import numpy as np
import json
import streamlit as st
import os
import csv
import logging
import re
import subprocess
import sys
import argparse
import time
import io
import shutil

st.set_page_config(page_title="Egcore study Validation and Convertion", page_icon="🧩")



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
    st.write(f"## Validation Based on CV version `{cv['version']}`.")
    ddsyn2={}
    for x in cv['classes'][1]['terms'].keys():
        try:
            temp=cv['classes'][1]['terms'][x]["alternative_terms"].strip().split(",")
            #print(temp)
            temp.append(cv['classes'][1]['terms'][x]["cv_name"])
            temp.append(cv['classes'][1]['terms'][x]["display_name"])
            #print(temp)
        except:
            temp=[]
            temp.append(cv['classes'][1]['terms'][x]["cv_name"])
            temp.append(cv['classes'][1]['terms'][x]["display_name"])
        
        try:
            temp2=cv['classes'][1]['terms'][x]["synonyms"].strip().split(",")
        except:
            temp2=[]
            
        temp.extend(temp2)
        
        temp=[y.strip() for y in temp if len(y) > 0]
        temp=[y.strip() for y in temp if len(y) > 0]
        for zz in temp:
            ddsyn2[zz] = cv['classes'][1]['terms'][x]["cv_name"]

    """# ### Load Entity and SIM mapping"""
    ddEntity={cv['classes'][1]['terms'][x]["cv_name"]:cv['classes'][1]['terms'][x]['parent_entity'] for x in cv['classes'][1]['terms'].keys()}
    ddSIM={cv['classes'][1]['terms'][x]["cv_name"]:cv['classes'][1]['terms'][x]['display_name'] for x in cv['classes'][1]['terms'].keys()}
    ddDType={cv['classes'][1]['terms'][x]["cv_name"]:cv['classes'][1]['terms'][x]['datatype'] for x in cv['classes'][1]['terms'].keys()}
    
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
                
    ddNumRange={cv['classes'][1]['terms'][x]["cv_name"]:eval(ddNumRange[x]) for x in ddNumRange}
    ddCatRange={cv['classes'][1]['terms'][x]["cv_name"]:[_.strip() for _ in ddCatRange[x].split(',')] for x in ddCatRange}
        
    return ddsyn2, ddSIM, ddEntity, ddDType, ddNumRange, ddCatRange

### Download data from S3 bucket
def createMappingURL(sid, outpath="./"):
    """Create mapping.csv file urls from study id and the output file path"""
    url = f"aws s3 cp s3://eaglecore-unilever/studies/{sid}/mapping.csv {outpath}{sid}_mapping.csv"
    return url

def createStudyDetailsURL(sid, outpath="./"):
    """Create StudyDetails.json file urls from study id and the output file path"""
    url = f"aws s3 cp s3://eaglecore-unilever/studies/{sid}/StudyDetails.json {outpath}{sid}_StudyDetails.json"
    return url

### Load mapping.csv file func
def LOAD_MAPPING_FILE(SID, ddsyn2, delim):
    """# ### Load mapping.csv file into a pandas DataFrame.
    process the DataFrame and map the column headers.
    Remove duplicate columns
    Drops Blank columns.

    return: a DataFrame
    """
    x=f"./egcoreInputs/egcore{SID}/{SID}_mapping.csv"    ##FILE LOCATION
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
    st.write(f"DataFrame size: {mappingFile.shape}")

    # ### Delete Blank columns
    mappingFile = mappingFile.dropna(axis=1, how='all')
    logging.info("DataFrame: BLANK Attributes removed.")
    logging.info(f"DataFrame new size: {mappingFile.shape}")
    
    # ### strip all values
    mappingFile = mappingFile.apply(lambda x: x.apply(lambda y: y.strip() if type(y) == type('') else y), axis=0)
    return mappingFile

def convertDtype(dseries, ddDType, errors="raise"):
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
        dseries = dseries.apply(lambda x: pd.to_datetime(x, format='%d_%m_%Y', errors=errors)).dt.date
        #logging.info(f"converted dataType of {s}: {dt2} --> DateTime, error mode = {errors}")
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
            #st.write(f"{s} have out of range numerical values")
            logging.error(f"Atribute {Ndata.name}: contains Out of Range Numerical Values")
            st.write(ddNumRange[s])
            st.write(Ndata)
            outOfRngCol.append(Ndata.name)
    if dt == 'String' and s in ddCatRange:
        Cdata=dseries.copy()
        if Cdata.dropna().isin(ddCatRange[s]).all():
            pass
        else:
            #st.write(f"{s} have out of range categorical values")
            logging.error(f"Atribute {Cdata.name}: contains Out of Range Categorical Values")
            st.write(ddCatRange[s])
            st.write(Cdata)
            outOfRngCol.append(Cdata.name)
    
    return outOfRngCol

def checkValidNumValues(dseries, ddNumRange, ddSIM):
    s=dseries.name
    outOfRngCol=[]
    
    #fp="/Users/raktimmaiti_mbp/Documents/GitHub/EagleGenomics-CV/master/eaglegenomics-cv.json"
    #with open(fp) as f:
        #cv=json.load(f)
    
    if s in ddNumRange:
        Ndata=dseries.copy()
        #cvname=cv['classes'][1]['terms'][f'{Ndata.name}_Study']["display_name"]
        cvname=ddSIM[Ndata.name]
        if Ndata.dropna().between(ddNumRange[s][0],ddNumRange[s][1]).all():
            pass
        else:
            #st.write(f"{s} have out of range numerical values")
            logging.error(f"Atribute {cvname}: contains Out of Range Numerical Values")
            st.write(cvname, "Valid Range: ", ddNumRange[s])
            st.write(cvname, "Actual data-> min: ", Ndata.min(),", max: ", Ndata.max())
            outOfRngCol.append(cvname)
    else:
        pass
    
    return outOfRngCol


def checkValidCatValues(dseries, ddCatRange, ddSIM):
    s=dseries.name
    outOfRngCol=[]
    
    #fp="/Users/raktimmaiti_mbp/Documents/GitHub/EagleGenomics-CV/master/eaglegenomics-cv.json"
    #with open(fp) as f:
        #cv=json.load(f)
    
    if s in ddCatRange:
        Cdata=dseries.copy()
        #cvname=cv['classes'][1]['terms'][f'{Cdata.name}_Study']["display_name"]
        cvname=ddSIM[Cdata.name]
        if Cdata.dropna().isin(ddCatRange[s]).all():
            pass
        else:
            #st.write(f"{s} have out of range categorical values")
            logging.error(f"Atribute {cvname}: contains Out of Range Categorical Values")
            st.write(cvname, "Valid Values: ", ddCatRange[s])
            #st.write(Cdata.name, "Actual Values: ", Cdata.unique())
            st.write(cvname, "Unknown Values: ", set(Cdata.unique())-set(ddCatRange[s]))
            outOfRngCol.append(cvname)
    else:
        pass
    
    return outOfRngCol


# ## Load BLANK Templates
def LOAD_BLANK_TP():
    """# ## Load BLANK Templates.
    returns a list of DFs
    """
    meta=pd.read_excel('./input_blank/20230412_Onboarding_Template_Default_V3_forBackendUse.xlsx', sheet_name=None, header=0)
    metacopy=meta.copy()
    #logging.info("BLANK MetaData TP loaded.")

    ########## DEFINE ALL SHEETS as DFs ##########
    study= meta['Study_Experiment']
    sample=meta['Data']

    return study, sample


def convert_df(df):
    return df.to_csv().encode('utf-8')



#########################################################################################################

        
        
        
def main():
    
    st.header("Eagle Core Study Validation and convertion!")
    
    loglevel=10 ## DEBUG level
    SID = st.text_input('EagleCore Study ID:', key="SID")
    errors = st.radio('Errors handled as:', ['raise', 'coerce', 'ignore'], horizontal=True, key="errors")
    #delim = st.text_input(label='The delimiter to be used:', value=",")
    delim = st.selectbox(label='The delimiter to be used:', options=[",", ";", "\t"], key="delim", help=None)
    
    submitted=0
    
    submit_button = st.button(label='Submit and Fetch files')
    #submit=st.select_slider(label='Submit and Fetch files', options=("No", "Yes"), 
    #                 key=None, help="Submit StudyID and Fetch files from Eaglecore S3 bucket")
    
    st.write("----------------------------------------------------------------------------------------------")
    st.write("Study ID: ", SID)
    st.write("Level of logging: ", loglevel)
    
    
    #timestr = time.strftime("%Y-%m-%d_%H_%M_%S")
    dir=f'./egcoreInputs/egcore{SID}/log'
    try:
        os.makedirs(dir)
        st.write(f"Folder {dir} created.")
    except Exception:
        st.write("Folder already exists.")
    logging.basicConfig(filename=f'./egcoreInputs/egcore{SID}/log/{SID}.log', 
                        format='%(levelname)s>> %(asctime)s: %(message)s', 
                        datefmt='%d/%m/%Y %I:%M:%S %p',
                        encoding='utf-8', level=loglevel)
    
    logging.info(f'STARTED at: {time.strftime("%a, %d %b %Y %H:%M:%S")}')
    buf = io.StringIO()
    
    #st.error('Error message')
    #st.warning('Warning message')
    #st.info('Info message')
    #st.success('Success message')
    
    
    
    if submit_button:
        submitted+=1
    
    # st.form_submit_button returns True upon form submit
    if submitted>0:
        with open("egcore_downloadfiles.sh","w") as f1:
            m=createMappingURL(SID, outpath=f"./egcoreInputs/egcore{SID}/")
            s=createStudyDetailsURL(SID, outpath=f"./egcoreInputs/egcore{SID}/")
            f1.write(m+'\n')
            f1.write(s)
        
        r=subprocess.call(['sh', './egcore_downloadfiles.sh'])
        if r==0:
            logging.info("Files from S3 bucket downloaded successfully.")
            st.success("Files from S3 bucket downloaded successfully.")
        else:
            logging.critical("Files from S3 bucket could not be downloaded.")
            st.error("Files from S3 bucket could not be downloaded.")
        ## Delete the temp file
        os.remove("./egcore_downloadfiles.sh")
    
    
    
    validate = st.checkbox('Validate and Curate')
    if validate:
        
        dataValidation=[]

        ddsyn2, ddSIM, ddEntity, ddDType, ddNumRange, ddCatRange = CREATE_VALIDATION_DICTS()
        mappingFile = LOAD_MAPPING_FILE(SID, ddsyn2, delim=delim)
        try:
            mappingFile = mappingFile.drop(['specimenDescription'], axis=1)
        except:
            ''
    
    
    
    
        study, sample = LOAD_BLANK_TP()
     
        ###STUDY EXPERIMENT FIXED INFO
        
        #try:
        SEfile="/Users/raktimmaiti_mbp/Library/CloudStorage/OneDrive-SharedLibraries-EagleGenomics/Services - Documents/Data_Curation_working_Folder/Eaglecore/EGCORE_STUDY_EXP_master/EGCORE_StudyName_Mapping_REVISED3.xlsx"
        SEInput=pd.read_excel(SEfile, sheet_name="Study_Experiment_REVISED_NEW", index_col="Study Identifier")
        SID=int(SID)
        StudyDict={}
        StudyDict["Study Identifier"]=SID
        StudyDict["Study Title"]=SEInput.loc[SID]['STUDY_TITLE']
        StudyDict["Study Type"]=SEInput.loc[SID]['studyType']
        StudyDict["ELN"]=SEInput.loc[SID]['eln']
        StudyDict["Experiment Name"]=SEInput.loc[SID]['experimentName']
        StudyDict["Experiment Technology"]=SEInput.loc[SID]['Assay Technology Type']
        StudyDict["Experiment Type"]=SEInput.loc[SID]['experimentType']
        StudyDict["Platform"]=SEInput.loc[SID]['platform']
        StudyDict["InvestigationType"]=SEInput.loc[SID]['investigationType']
        StudyDict["PairedOrSingle"]=SEInput.loc[SID]['Indexed Sequencing']
        StudyDict["Processing Status"]=SEInput.loc[SID]['processingStatus']
        StudyDict["Raw Data Available"]=SEInput.loc[SID]['rawDataAvailable']
        
        #except:
        #    StudyDict={}
        #    StudyDict["Study Identifier"]=SID
        #    StudyDict["Study Type"]="TEST_Oral"
        #    StudyDict["ELN"]="TEST_EXP-17-AH2889"
        #    StudyDict["Experiment Name"]="TestGAPFILLER"
        #    StudyDict["Experiment Technology"]="NGS"
        #    StudyDict["Experiment Type"]="Metagenomics"
        #    StudyDict["Platform"]="Illumina"
        #    StudyDict["InvestigationType"]="16S"
        #    StudyDict["Processing Status"]="Raw"
        #    StudyDict["Raw Data Available"]="Yes"
        
        
        
        
        #################################################################
        
        SSstudy=pd.DataFrame.from_dict(StudyDict, orient='index').T
        
        
        allDefined=[x for x in mappingFile.columns if "~" in x or "_|" in x]
        if len(allDefined) > 0:
            logging.critical("mapping file have Undefined or Duplicate attributes.")
            st.error(f"Undefined or Duplicate attributes present: {allDefined}")
            st.write(mappingFile.head(5))
            #sys.exit("mapping.csv file have Undefined or Duplicate attributes.")
            dataValidation.append(False)
        else:
            dataValidation.append(True)
            
            
        #fp="/Users/raktimmaiti_mbp/Documents/GitHub/EagleGenomics-CV/master/eaglegenomics-cv.json"
        #with open(fp) as f:
        #    cv=json.load(f)
            
        
        
        SE=pd.DataFrame()
        DATA=pd.DataFrame()
        
        OutOfRngCols=[]
        InvalidDataTypeCols=[]
        
        for columns in mappingFile:
            x=mappingFile[columns]
            try:
                x= convertDtype(x, ddDType, errors=errors)
            except ValueError:
                InvalidDataTypeCols.append(ddSIM[x.name])
                continue
            
            #valChk=checkValidValues(x, ddDType, ddNumRange, ddCatRange)
            #OutOfRngCols.extend(valChk)
            
            if columns in ddSIM:
                sim = ddSIM[columns]
                entryPoint=ddEntity[columns]
                dataType=ddDType[columns]
            else:
                sim = columns
                entryPoint = 'Sample'
                dataType= 'String'
                
            if dataType == "Numeric" or dataType == "INT":
                valChk=checkValidNumValues(x, ddNumRange, ddSIM)
            elif dataType == "String":
                valChk=checkValidCatValues(x, ddCatRange, ddSIM)
            
            OutOfRngCols.extend(valChk)
            
            
            prefix=[]
            data = [*prefix, *x]
            try:
                dname=ddSIM[columns]
            except:
                dname=columns
            data = pd.Series(data, name=dname)
        
            #st.write(data)
            
            #### Fill other tabs
            if entryPoint == 'Experiment':
                SE = pd.concat([SE,data],axis=1)
        
            elif entryPoint == 'Study':
                SE = pd.concat([SE,data],axis=1)
            else:
                DATA = pd.concat([DATA,data],axis=1)
       
    
        # ### Check Out of Range Value columns
        if len(OutOfRngCols) > 0:
            st.error(f"{str(OutOfRngCols)} attributes have out of range values.")
            logging.critical(f"{str(OutOfRngCols)} attributes have out of range values.")
            dataValidation.append(False)
        else:
            dataValidation.append(True)
            logging.info("All mapped attributes have values within range wherever applicable.")
            
        # ### Check Invalid Value columns
        if len(InvalidDataTypeCols) > 0:
            st.error(f"{str(InvalidDataTypeCols)} attributes have Invalid Type of Data.")
            logging.critical(f"{str(InvalidDataTypeCols)} attributes have out of range values.")
            dataValidation.append(False)
            #sys.exit(f"{str(InvalidDataTypeCols)} attributes have Invalid Data.")
        else:
            logging.info("All mapped attributes have values with Valid Type of Data.")
            dataValidation.append(True)        
        
        
        DATA.insert(1, "Organism", ["Human"]*(len(DATA['SampleID'])))
        
        SE=SE.drop_duplicates()
        
        # ### Merge study tabs
        study2 = pd.concat([SE,SSstudy], axis=1)
        study2 = study2.loc[:,~study2.columns.duplicated()].copy()
        study2 = study2.dropna(axis=0, how='all')
        
        cols = study.columns.to_list()
        study2 = study2[cols + [c for c in study2.columns if c not in cols]]
        if len(study2["Study Identifier"]) >1:
            st.warning("WARNING: MULTIPLE ENTRY IN STUDY")
            st.write(study2)
            logging.critical("MULTIPLE ROWS of entry found in study tab.")
        
    
        
        
        
        ### Check all validation
        validationDict=dict(zip(["No New/Syn Attributes", 
                                 #"Mandatory Attributes Present", 
                                 "Valid Range", 
                                 "Valid Data Type"], dataValidation))
        if all(dataValidation):
            logging.info(f"ALL DATA VALIDATION PASSED. Details: {validationDict}")
            st.success(f"ALL DATA VALIDATION PASSED.")
            st.write(validationDict)
        else:
            logging.error(f"DATA VALIDATION FAILED. Details: {validationDict}")
            st.error(f"DATA VALIDATION FAILED.")
            st.write(validationDict)
        
        
        #col1, col2 = st.columns(2)
        
        st.caption("Study_Experiment")
        st.dataframe(study2)
        st.caption("Data")
        st.dataframe(DATA)
    

        study2 = study2.columns.to_frame().T.append(study2, ignore_index=True)
        study2.columns = range(len(study2.columns))
        DATA = DATA.columns.to_frame().T.append(DATA, ignore_index=True)
        DATA.columns = range(len(DATA.columns))
        
        
        dir=f'./CURATED/{SID}'
        try:
            os.makedirs(dir)
            st.write(f"Folder {dir} created.")
        except Exception:
            st.write("Folder already exists.")
        
        outfile=f"./CURATED/{SID}/EGCORE_Unilever_{SID}_CuratedMetaData.xlsx"
        
        with pd.ExcelWriter(outfile) as writer:
            study2.to_excel(writer, sheet_name='Study_Experiment', header=True, index=False)
            DATA.to_excel(writer, sheet_name='Data', header=True, index=False)
        
        try:
            with open(outfile, 'rb') as f:
                st.download_button('Download Curated Metadata File',
                                   f, file_name=f'EGCORE_Unilever_{SID}_CuratedMetaData.xlsx')
        except:
            ''
        
        #try:
        #    with open(f'./egcoreInputs/egcore{SID}/log/{SID}.log', 'r') as f2:
        #        logf=f2.read()
        #        st.code(logf, language="log", line_numbers=False)
        #except:
        #    ''
    
if __name__ == '__main__':
    main()
