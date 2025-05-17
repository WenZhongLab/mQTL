# Risk group prediction


from xgboost import XGBClassifier
from sklearn.metrics import roc_curve, auc
import pandas as pd
import numpy as np  

# use a variable list for predictor selection
def xbgROC_risk(name,refn,plfn):
    import warnings
    warnings.filterwarnings("ignore")
    roclist = []
    roclist.append( ".".join(name))
    cc=[]
    truelist = []
    problist = []
    fwb = open(refn,'w') 
    fwplot = open(plfn, 'w')
    for i in range(1,101):

        # load training set
        fn = "".join(['riskscale_train',"_time_",str(i),".csv"])
        train = pd.read_csv(fn)
        x_train = train.loc[:,name]
        y_train = train.loc[:,["group"]]
        # load validation set
        fn = "".join(['riskscale_test',"_time_",str(i),".csv"])
        test = pd.read_csv(fn)
        x_test = test.loc[:,name]
        y_test = test.loc[:,["group"]]
    
        seedx=i
        xgb = XGBClassifier(
                    objective = 'binary:logistic',
                    use_label_encoder=False,
                    subsample=0.5,  
                    scale_pos_weight=235/113,
                    random_state=seedx)
        xgb.fit(x_train, y_train, eval_metric = "auc")
        y_prob = xgb.predict_proba(x_test)[:,1]
        fpr, tpr, thresholds = roc_curve(y_test, y_prob, pos_label=1)
        c = auc(fpr,tpr)
        roclist.append(str(format(c, '.6f')))
        cc.append(c)
        truelist.extend(np.transpose(y_test).values.tolist()[0])
        problist.extend(y_prob.tolist())

    # export auc
    fwb.write(",".join(roclist))
    fwb.write("\n")
    
    # export true group and probality for ROC
    truelistnew = [str(x) for x in truelist]
    problistnew = [str(x) for x in problist]
    fwplot.write(",".join(truelistnew))
    fwplot.write("\n")
    fwplot.write(",".join(problistnew))
    fwplot.write("\n")

    if(i==100):
        print(name)
        print("Finish 100 times prediction")
        print(np.mean(cc))
    fwb.close()
    fwplot.close()


# use index of the start and end of predicters in the data set
def xbgROC_risk_index(name,name_start,name_end,refn, plfn):
    import warnings
    warnings.filterwarnings("ignore")

    roclist = []
    roclist.append( ".".join(name))
    cc=[]
    truelist = []
    problist = []
    fwb = open(refn,'w') 
    fwplot = open(plfn, 'w')
    for i in range(1,101):
        # load training set
        fn = "".join(['riskscale_train',"_time_",str(i),".csv"])
        train = pd.read_csv(fn)
        x_train = train.loc[:,name_start:name_end]
        y_train = train.loc[:,["group"]]
        # load validation set
        fn = "".join(['riskscale_test',"_time_",str(i),".csv"])
        test = pd.read_csv(fn)
        x_test = test.loc[:,name_start:name_end]
        y_test = test.loc[:,["group"]]

        seedx=i

        xgb = XGBClassifier(
                    objective = 'binary:logistic',
                    use_label_encoder=False,
                    subsample=0.5,  
                    scale_pos_weight=235/113,
                    random_state=seedx)
        xgb.fit(x_train, y_train, eval_metric = "auc")
        y_prob = xgb.predict_proba(x_test)[:,1]
        fpr, tpr, thresholds = roc_curve(y_test, y_prob, pos_label=1)
        c = auc(fpr,tpr)
        roclist.append(str(format(c, '.6f')))
        cc.append(c)
        truelist.extend(np.transpose(y_test).values.tolist()[0])
        problist.extend(y_prob.tolist())


    # export auc
    fwb.write(",".join(roclist))
    fwb.write("\n")
    
    # export true group and probality for ROC
    truelistnew = [str(x) for x in truelist]
    problistnew = [str(x) for x in problist]
    fwplot.write(",".join(truelistnew))
    fwplot.write("\n")
    fwplot.write(",".join(problistnew))
    fwplot.write("\n")

    if(i==100):
        print(name)
        print("Finish 100 times prediction")
        print(np.mean(cc))
    fwb.close()
    fwplot.close()


# run function, use a variable list for predictor selection
name=["ANGPTL4","LEP","TMPRSS15","FAM3B","FAM3C","LDLR","FGF21","DCBLD2","IL1RN","CTSF","SERPINA9", "CPM" ,"GCG"]
xbgROC_risk(name,"./result/tier1protein_100balance_pro13.csv",
               "./result/tier1protein_100balance_pro13scale_forROC.csv")


name=["lc_132", "lc_122", "lc_115", "gc_47" , "gc_94",  "gc_166", "gc_86" , "lc_158" ,"gc_65" , "lc_185" ,"gc_74",  "lc_178", "lc_244" ,"lc_300"]
xbgROC_risk(name,"./result/tier1met_100balance_met14.csv",
               "./result/tier1met_100balance_met14scale_forROC.csv",)


name=["ANGPTL4","LEP","TMPRSS15","FAM3B","FAM3C","LDLR","FGF21","DCBLD2","IL1RN","CTSF","SERPINA9", "CPM" ,"GCG","lc_132", "lc_122", "lc_115", "gc_47" , "gc_94",  "gc_166", "gc_86" , "lc_158" ,"gc_65" , "lc_185" ,"gc_74",  "lc_178", "lc_244" ,"lc_300"]
xbgROC_risk(name,"./result/tier1promet_100balance_pro13met14.csv",
               "./result/tier1promet_100balance_pro13met14scale_forROC.csv")


# run function, use index of the start and end of predicters in the data set
xbgROC_risk_index("allprotein","ACAN","ZBTB17",
                     "./result/allprotein_100balance_pro794.csv",
                     "./result/allprotein_100balance_pro794_forROC.csv")


xbgROC_risk_index("allmet","gc_18","lc_272",
                     "./result/allmet_100balance_met527.csv",
                     "./result/allmet_100balance_met527_forROC.csv")


xbgROC_risk_index("allpromet","ACAN","lc_272",
                     "./result/allproteinmet_100balance_pro794met527scale.csv",
                     "./result/allproteinmet_100balance_pro794met527_forROC.csv")




