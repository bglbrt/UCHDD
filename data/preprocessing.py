# LIBRARIES
import pandas as pd
import numpy as np
import os

DATA_PATH = "PATH/TO/DATA/FOLDER/"
RAW_DATA_PATH = "PATH/TO/RAW/DATA/FOLDER/"

################################################################################
# PART III ESSAY - UNMEASURED CONFOUNDING IN HIGH-DIMENSIONAL DATA
# - SIMULATIONS
# - FILE N.7BIS
################################################################################

################################################################################
# PREPROCESSING

# list all TSV files (files from main study)
files = []
for file in os.listdir(RAW_DATA_PATH):
    if file.endswith(".tsv"):
        files = files + [file]

# append all covariates in a single dataset
df = pd.read_csv(os.path.join(RAW_DATA_PATH, files[0]), sep='\t')
new_columns = [str(col).strip() for col in df.columns.values]
new_columns[0] = "temp"
df.columns = new_columns
df = df[["temp", "2011"]]
df[['var_code', 'CITY']] = df["temp"].str.split(',', 1, expand=True)
df = df[["CITY", "var_code", "2011"]]
df = df.loc[df["CITY"].str[0:2] == "DE"]
df.head(5)
df = df.pivot(index=["CITY"], columns=["var_code"])
df.columns = df.columns.droplevel(0)
df.columns.name = None
df = df.reset_index()
df = df[df["CITY"] != "DE"]
for i in range(1, len(files)):
    _ = pd.read_csv(os.path.join(RAW_DATA_PATH, files[i]), sep='\t')
    new_columns = [str(col).strip() for col in _.columns.values]
    new_columns[0] = "temp"
    _.columns = new_columns
    _ = _[["temp", "2011"]]
    _[['var_code', 'CITY']] = _["temp"].str.split(',', 1, expand=True)
    _ = _[["CITY", "var_code", "2011"]]
    _ = _.loc[_["CITY"].str[0:2] == "DE"]
    _.head(5)
    _ = _.pivot(index=["CITY"], columns=["var_code"])
    _.columns = _.columns.droplevel(0)
    _.columns.name = None
    _ = _.reset_index()
    _ = _[_["CITY"] != "DE"]
    df = df.merge(_, how='inner', left_on="CITY", right_on="CITY")

# add environment variables from 2012 study (posits that this has not changed)
# covariates include primarily the shares of land area in terms of:
# - transport infrastructure
# - residential urban areas
# - agricultural areas
# - natural areas
# - industrial areas
_ = pd.read_csv(os.path.join(RAW_DATA_PATH, "urb_cenv.tsv"), sep='\t')
new_columns = [str(col).strip() for col in _.columns.values]
new_columns[0] = "temp"
_.columns = new_columns
_ = _[["temp", "2012"]]
_[['var_code', 'CITY']] = _["temp"].str.split(',', 1, expand=True)
_ = _[["CITY", "var_code", "2012"]]
_ = _.loc[_["CITY"].str[0:2] == "DE"]
_.head(5)
_ = _.pivot(index=["CITY"], columns=["var_code"])
_.columns = _.columns.droplevel(0)
_.columns.name = None
_ = _.reset_index()
_ = _[_["CITY"] != "DE"]
_ = _[[col for col in _.columns if col[0:4] in ["EN52", "CITY"]]]

# replace all different NaN value with a single NaN value
df = df.replace(': ', np.nan)
df = df.replace(':', np.nan)
df = df.replace(':', np.nan)
df = df.replace('NaN', np.nan)
df = df.replace('NaN ', np.nan)
df = df.replace('nan', np.nan)

# drop cities with NaN in all columns
df = df.dropna(axis=1, how='all')
df = df.merge(_, how='inner', left_on="CITY", right_on="CITY")

# convert all columns to float (some columns were integers although continuous variables)
for col in df.columns[1:]:
    df[col] = [float(n) if n not in ["nan", ":", ": ", "NaN", "nan"] else np.nan for n in df[col].astype(str).str.split(' ', expand=True)[0]]

# add cities names and geographical information (latitude and longitude)
df = df[df["CR2001V"].isna() == False]
codes = pd.read_csv(os.path.join(RAW_DATA_PATH, "urb_esms.csv"), sep=',')
df = df.merge(codes, how="inner", left_on="CITY", right_on="CODE")
df = df[["NAME"] + sorted(list(df.columns[1:-2]))]
latlong = pd.read_csv(os.path.join(RAW_DATA_PATH, "urb_geo.csv"), sep=',')
df = df.merge(latlong, how="inner", left_on="NAME", right_on="CITY")
df = df[["NAME", "LAT", "LONG"] + sorted(list(df.columns[1:-3]))]

# set NaN to column mean for the (few) missing data
column_means = df.mean()
df = df.fillna(column_means)

# export dataset
df.to_csv(DATA_PATH+"data.csv")
