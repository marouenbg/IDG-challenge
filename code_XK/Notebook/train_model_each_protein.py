import math
import random
import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold
from paths import path_out_data, path_data
import GPy
from scipy.stats import pearsonr,spearmanr
from sklearn.metrics import mean_squared_error
from sklearn.linear_model import ElasticNet
import sys
from glob import glob
from tqdm import tqdm
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV
from sklearn.kernel_ridge import KernelRidge
from sklearn.svm import SVR
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.metrics import mean_squared_error as rmse
from sklearn.ensemble import RandomForestRegressor

index_feature = int(sys.argv[1])
index_model = int(sys.argv[2])
version=sys.argv[3]
sanity = False


train_test = pd.read_csv("%s/train_test/SMILE_prot_id_kd.csv"%(path_out_data))
file_features = glob("%s/unified_features_compound/*_ki.csv"%(path_out_data))
print("feature set in  total: ", len(file_features))
print("Invalid Kd", sum(train_test.Kd < 0 ))
#you might see warnings becuase the testing rows don't have Kd
random.Random(4).shuffle(file_features)
print(file_features)
file_feature = file_features[index_feature]
feature = pd.read_csv(file_feature)
train_test = pd.merge(train_test,feature, how='left', on="SMILE")
print(file_feature)

mask = train_test.Kd.isna()
test = train_test[mask]
protein_u_t = test.target_id.unique() #u:unique, t:test


num_tr = []
num_te = []
protein_missing = []
for protein in protein_u_t:
    try: 
        data = train_test[train_test.target_id == protein]
        data = data.loc[:,data.apply(pd.Series.nunique) != 1]
        mask = data.Kd.isna()
        train = data[~mask]
        test = data[mask]
        num_tr.append(train.shape[0])
        num_te.append(test.shape[0])
        if train.shape[0] ==0:
            protein_missing.append(protein)
    except:
        print(protein)
        protein_missing.append(protein)
        
protein_u_t = [p for p in protein_u_t if p not in protein_missing]


def train_SVR(X_fold_tr,Y_fold_tr):
    svr = GridSearchCV(SVR(kernel='rbf'), cv=3,
                       param_grid={"C": [1e0, 1e1, 1e2, 1e3],
                               "gamma": np.logspace(-2, 2, 5)},
                      scoring='neg_mean_absolute_error', n_jobs=-1)
    svr.fit(X_fold_tr,Y_fold_tr)
    model = svr.best_estimator_  
    model.fit(X_fold_tr, Y_fold_tr)
    
    return model

def train_kernel_ridge(X_fold_tr,Y_fold_tr):
    """
    run cv on the training fold and return the best model
    """
    
    KRmodel = GridSearchCV(KernelRidge(), cv = 3,
              param_grid={"alpha": np.logspace(-10, -5, 10),
             "gamma": np.logspace(-12, -9, 10), 
             "kernel" : ['laplacian', 'rbf']}, 
             scoring='neg_mean_absolute_error', n_jobs=-1)

    KRmodel = KRmodel.fit(X_fold_tr, Y_fold_tr)
    model = KRmodel.best_estimator_        
    model.fit(X_fold_tr, Y_fold_tr) 
    return model

def train_RF(X_fold_tr, Y_fold_tr):
    
    RFmodel = GridSearchCV(RandomForestRegressor(), cv= 3,
              param_grid={"n_estimators": np.linspace(50, 150, 25).astype('int')}, 
                           scoring='neg_mean_absolute_error', n_jobs=-1)

    RFmodel = RFmodel.fit(X_fold_tr, Y_fold_tr)
    model = RFmodel.best_estimator_        
    model.fit(X_fold_tr, Y_fold_tr) 
    return model


def test_model(X_train, Y_train,X_test,model_fun):
    """
    Use kernel ridge regression
    
    """
    
    preds = []
    truths = []
    indexes = []
    pred_test = np.zeros(X_test.shape[0])    
    
    bins = np.linspace(0,  max(Y_train),3)
    y_binned = np.digitize(Y_train, bins)    
    num_fold = 3
    skf = StratifiedKFold(n_splits=num_fold)

    
    for train_index, valid_index in skf.split(X_train, y_binned):

        X_fold_tr = X_train.iloc[train_index].values
        #mask for constant zero if any
        mask = ~(np.apply_along_axis(lambda x: max(x)-min(x), 0, X_fold_tr) ==0)
        X_fold_tr = X_train.loc[:,mask].iloc[train_index]
        Y_fold_tr = Y_train.iloc[train_index].values
        X_fold_valid = X_train.loc[:,mask].iloc[valid_index].values
        Y_fold_valid = Y_train.iloc[valid_index].values
        
        #run CV on the fold for training
        model = model_fun(X_fold_tr,Y_fold_tr)
        
        # use the best model to predict on the fold_valid
        pred  = model.predict(X_fold_valid).flatten()
        pred_test += model.predict(X_test.loc[:,mask]).flatten()

        preds.append(pred)
        truths.append(Y_fold_valid)
        indexes.append(X_train.iloc[valid_index].index)
    pred_test /=num_fold
    
    truths = np.concatenate(truths,axis=0).flatten()
    preds = np.concatenate(preds,axis =0 ).flatten()
    indexes = np.concatenate(indexes).flatten()
    return(preds,truths, indexes, pred_test)


model_funs = [train_SVR, train_kernel_ridge, train_RF ]
model_fun = model_funs[index_model]


pred_all_proteins = []
truth_all_proteins =[]
index_all_proteins =[]
SMILE_protein_tests = []

if sanity:
    pbar = tqdm(enumerate(protein_u_t[:10]), total = len(protein_u_t[:10]))
else:
    pbar = tqdm(enumerate(protein_u_t), total = len(protein_u_t))

for i,protein in pbar:
    data = train_test[train_test.target_id == protein]
    data = data.loc[:,data.apply(pd.Series.nunique) != 1]  
    
    features = data.loc[:, data.columns != 'Kd']
    # set(features.dtypes) 
    # exclude non-numerical columns
    numerics = ['int64','float64']
    features = features.select_dtypes(include=numerics)
    #normalize features
    features = features.apply(lambda x: (x-x.min())/(x.max()-x.min()), axis=0)
       
    if features.shape[0] < features.shape[1]:
        try:
            pca = PCA(n_components= min(features.shape))
            features = pd.DataFrame(pca.fit_transform(features.values))  
            features.index = data.index
        except:
            pass

    mask = data.Kd.isna()
    X_train = features[~mask]
    Y_train = data[~mask].Kd
    X_test = features[mask]
    # if the number training data is less than 3. skf.split will spit an error as it 3-fold
    if (X_train.shape[0]<3):
        continue
    try:
        preds, truths, indexes, pred_test = test_model(X_train, Y_train,X_test,model_fun)
    
        SMILE_protein_test = data[mask].iloc[:,:1]
        SMILE_protein_test["protein"] = protein
        SMILE_protein_test["pred"] = pred_test
    
        pred_all_proteins.append(preds)
        truth_all_proteins.append(truths)
        index_all_proteins.append(indexes)
        SMILE_protein_tests.append(SMILE_protein_test)
    except:
        print("too few samples for this protein %s"%(protein))
    
    

output = pd.DataFrame([np.concatenate(truth_all_proteins),
                       np.concatenate(pred_all_proteins),
                       np.concatenate(index_all_proteins)])
output = output.transpose()
output.columns = ["truth","preds","indexes"]
output_name = file_feature.split("/")[-1][:-4] + model_fun.__name__
output.to_csv("%s/predictions/%s_%s_CV.csv"%(path_out_data,output_name,version),index = False)

pred_test = pd.concat(SMILE_protein_tests,axis=0)
pred_test.to_csv("%s/predictions/%s_%s_test.csv"%(path_out_data,output_name,version),index = False)


PCC = pearsonr(np.concatenate(truth_all_proteins),np.concatenate(pred_all_proteins))
SRC = spearmanr(np.concatenate(truth_all_proteins),np.concatenate(pred_all_proteins))
print("PCC: %s"%(PCC[0]))
print("SRC: %s"%(SRC[0]))

error = rmse(np.concatenate(truth_all_proteins),np.concatenate(pred_all_proteins))


print("rmse: %s"%(error))


perf = pd.DataFrame([PCC[0], SRC[0], error])
perf.to_csv("%s/predictions/%s_%s_perf.csv"%(path_out_data,output_name,version),index = False)

assert sum([i==j for i,j in zip(train_test.loc[output.indexes,:].Kd, output.truth)]) == output.shape[0]
