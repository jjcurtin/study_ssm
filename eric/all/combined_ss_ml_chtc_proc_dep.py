

sns.set_style('ticks')
sns.set_context('paper',font_scale=2.5)
fig,axs = plt.subplots(1,4,figsize=(20,7))
for i,delta in enumerate(window_delta):
    #sns.pointplot(data=eval_df.query("window==@delta"),x='horizon',y='auc',hue='fit_type',ax=axs[i],native_scale=True)
    sns.pointplot(data=eval_df.query("window==@delta"),x='day',y='ba',hue='fit_type',hue_order=hue_ord,ax=axs[i],native_scale=True)
    axs[i].set_title(pane_titles[i])
    axs[i].set_xlim(-5,95)
    axs[i].set_ylim(.7,1.0)
    axs[i].set_yticks(np.arange(0.7,1.0,0.05))
    axs[i].set_xticks(np.arange(0,105,15))
    axs[i].set_xlabel('Day')
    if i==0:
        axs[i].set_ylabel('Balanced Accuracy')
    else:
        axs[i].set_ylabel('')
    axs[i].legend(title='Model type')
#fig.suptitle('Balanced accuracy by time step for different prediction windows')
fig.tight_layout()
#fig.savefig("plots/ss_ml_ba.pdf",bbox_inches='tight',facecolor='w')


sns.set_style('whitegrid')
sns.set_context('paper',font_scale=2)
fig = plt.figure(figsize=(20,7))
subfigs = fig.subfigures(1,4,wspace=0.05)

for i,delta in enumerate(window_delta):
    ax = subfigs[i].subplots()
    #sns.pointplot(data=eval_df.query("window==@delta"),x='horizon',y='aucpr',hue='fit_type',ax=axs[i],native_scale=True)
    sns.lineplot(data=eval_df.query("fit_type in @hue_ord & window==@delta"),x='day',y='aucpr',hue='fit_type',hue_order=hue_ord,ax=ax,linewidth=2)#,native_scale=True)
    ax.set_title(pane_titles[i])
    ax.set_xlim(-5,95)
    ax.set_ylim(.1,1)
    ax.set_yticks(np.arange(0,1.1,0.1))
    ax.set_xticks(np.arange(0,105,15))
    ax.set_xlabel('Day')
    if i==0:
        ax.set_ylabel('AUCPR')
    else:
        ax.set_ylabel('')
    if i==3:
        ax.legend(title='Model type')
    else:
        ax.get_legend().remove()

handles, labels = ax.get_legend_handles_labels()

plt.legend(handles=handles,ncol=4,loc='upper center',bbox_to_anchor=(-1.55,-0.18))
#fig.suptitle('AUCPR by time step for different prediction windows')
#fig.tight_layout()
#fig.savefig("plots/ss_ml_aucpr_frozen.png",bbox_inches='tight',facecolor='w')

sns.set_style('whitegrid')
sns.set_context('paper',font_scale=2.5)
hue_ord=['MAP','MLE','LR','XGB']
fig = plt.figure(figsize=(20,7))
subfigs = fig.subfigures(1,4,wspace=0.05)

for i,delta in enumerate(window_delta):
    ax = subfigs[i].subplots()
    #sns.pointplot(data=eval_df.query("window==@delta"),x='horizon',y='aucpr',hue='fit_type',ax=axs[i],native_scale=True)
    sns.lineplot(data=eval_df.query("fit_type in @hue_ord & window==@delta"),x='day',y='auc',hue='fit_type',hue_order=hue_ord,ax=ax,linewidth=2)#,native_scale=True)
    ax.set_title(pane_titles[i])
    ax.set_xlim(-5,95)
    ax.set_ylim(.45,1.)
    ax.set_yticks(np.arange(0.45,1.05,0.05))
    ax.set_xticks(np.arange(0,105,15))
    ax.set_xlabel('Day')
    if i==0:
        ax.set_ylabel('AUC')
    else:
        ax.set_ylabel('')
    if i==3:
        ax.legend(title='Model type')
    else:
        ax.get_legend().remove()

handles, labels = ax.get_legend_handles_labels()

plt.legend(handles=handles,ncol=4,loc='upper center',bbox_to_anchor=(-1.55,-0.18))
#fig.suptitle('AUCPR by time step for different prediction windows')
#fig.tight_layout()
#fig.savefig("plots/ss_ml_auc_frozen.png",bbox_inches='tight',facecolor='w')

for fit_type in ['MAP_mle','MLE']:
    if fit_type=='MAP_mle':
        title = 'MAP'
    else:
        title='MLE'
    confusion_df = df.query("train_horizon>=2 & fit_type==@fit_type")
    pred = confusion_df.w0_pred.dropna().to_numpy()
    act = confusion_df.w0_act.dropna().to_numpy()
    fpr,tpr,thresh = roc_curve(act,pred)
    opt = np.argmax(tpr-fpr)
    #print(opt)
    cut = thresh[opt]
    pred_thresh = pred>=cut
    disp=ConfusionMatrixDisplay(confusion_matrix(act,pred_thresh,normalize='true'))
    disp.plot()
    plt.title(title)
    plt.show()

for fit_type in ['MAP_mle','MLE']:
    if fit_type=='MAP_mle':
        title = 'MAP'
    else:
        title='MLE'
    confusion_df = df.query("train_horizon>=2 & fit_type==@fit_type")
    pred = confusion_df.w7_pred.dropna().to_numpy()
    act = confusion_df.w7_act.dropna().to_numpy()
    fpr,tpr,thresh = roc_curve(act,pred)
    opt = np.argmax(tpr-fpr)
    #print(opt)
    cut = thresh[opt]
    pred_thresh = pred>=cut
    disp=ConfusionMatrixDisplay(confusion_matrix(act,pred_thresh,normalize='true'))
    disp.plot()
    plt.title(title)
    plt.show()