# UK Biobank metabolic disease prediction

from xgboost import XGBClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc
import pandas as pd
import numpy as np  


def xbgROC_balance(disease,refn,plfn):
    import warnings
    warnings.filterwarnings("ignore")
    name = ["ANGPTL4","LEP","TMPRSS15","FAM3B","FAM3C","LDLR","FGF21","DCBLD2","IL1RN","CTSF","SERPINA9", "CPM" ,"GCG"]
    #li = name
    roclist = []
    roclist.append( ".".join(name))
    cc=[]
    truelist = []
    problist = []
    fwb = open(refn,'w')
    fwplot = open(plfn,'w')
    for i in range(1,101):
        # load balanced dataset for modeling (randomly selected controls with equal size of disease group at each time)
        fn = "".join([disease,"_time_",str(i),".csv"])
        data = pd.read_csv(fn)
        x = data.loc[:,name]
        y = data.loc[:,["group"]]

        seedx=i
        x_train, x_test,y_train, y_test = train_test_split(x, y,
                                                   test_size = 0.3,
                                                   stratify = y,
                                                   random_state = seedx)

        xgb = XGBClassifier(
                    objective = 'binary:logistic',
                    use_label_encoder=False,
                    subsample=0.5,
                    scale_pos_weight=1,
                    random_state=seedx)
        xgb.fit(x_train, y_train, eval_metric =   "auc")
        y_pred = xgb.predict(x_test).ravel()
        y_prob = xgb.predict_proba(x_test)[:,1]
        fpr, tpr, thresholds = roc_curve(y_test, y_prob, pos_label=1)
        c = auc(fpr,tpr)
        roclist.append(str(format(c, '.6f')))
        cc.append(c)
        truelist.extend(np.transpose(y_test).values.tolist()[0])
        problist.extend(y_prob.tolist())
    # export AUC
    fwb.write(",".join(roclist))
    fwb.write("\n")
    # export true group and probability for ROC
    truelistnew = [str(x) for x in truelist]
    problistnew = [str(x) for x in problist]
    fwplot.write(",".join(truelistnew))
    fwplot.write("\n")
    fwplot.write(",".join(problistnew))
    fwplot.write("\n")
    
    if(i==100):
        print(disease)
        print("Finish 100 times predictions")
        print(np.mean(cc))
    fwb.close()
    fwplot.close()


# run function
xbgROC_balance("T2D","T2D_100balance_pro13.csv","T2D_100balance_pro13_forAUC.csv")



