import numpy as np
import matplotlib.pyplot as plt
import matplotlib

def plotsdefault(xlab='',ylab='',\
                 lw=1,lfs=12,tfs=12,dpi=300,size_x=5,size_y=4,Grid=False):
    #plt.rcParams.update(plt.rcParamsDefault)
    # We like big plots, with big fonts
    plt.rcParams['axes.linewidth'] = lw
    #plt.rc('text', usetex=True)
    matplotlib.rcParams["figure.dpi"] = dpi
    #matplotlib.rc('font', family='sans-serif') 
    #matplotlib.rc('font', serif='Helvetica Neue') 

    fig = plt.figure(figsize=(size_x,size_y))
    ax = fig.add_subplot(111)

    ax.set_xlabel(xlab,fontsize=lfs)
    ax.set_ylabel(ylab,fontsize=lfs)

#     ax.tick_params(which='major',direction='in',width=1,length=13,right=True,top=True,pad=7)
#     ax.tick_params(which='minor',direction='in',width=1,length=10,right=True,top=True)
    if Grid:
        ax.grid()
    return fig,ax

def plot_lobster(gn, gn3s, gn3sH, giIO, giIO3s, giIO3sH, colors = ['teal', 'royalblue']):
    
    fig, ax = plotsdefault(lfs=13,tfs=13)
    
    ax.fill_between(giIO3s[:,0],giIO3s[:,1],giIO3sH[:,1],alpha=0.5,color=colors[0],zorder=2)
    ax.fill(giIO3sH[:,0],giIO3sH[:,1],alpha=0.5,color=colors[0],zorder=2)
    ax.fill_between(giIO[:,0],giIO[:,1],color=colors[0], alpha=1, label=r'IO',zorder=3) #limegreen


    ax.fill_between(gn3s[:,0],gn3s[:,1],gn3sH[:,1],alpha = 0.5,color=colors[1],zorder=0)
    ax.fill_between(gn[:,0],gn[:,1],alpha = 1,color=colors[1],label=r'NO',zorder=1) #cornflowerblue
    
    effective_masses = np.ones(len(giIO3s[:,0]))
    lightest_neutrino = np.linspace(1e-4,1,len(giIO3s[:,0]))

    ax.text(9e-2, 6e-2, 'IO', color='white', fontsize=14, rotation=0, zorder=5)
    ax.text(1.3e-4, 2.0e-3, 'NO', color='white', fontsize=14, rotation=0, zorder=5)
   #ax.fill_between(lightest_neutrino,61e-3*effective_masses,165e-3*effective_masses,  color="silver", alpha=0.3,zorder=4)
   #ax.fill_between(lightest_neutrino,36e-3*effective_masses,165e-3*effective_masses,  color="silver", alpha=0.3,zorder=4) ## New limit 2022 KamLand Zen https://arxiv.org/abs/2203.02139
    #ax.text(1.5e-3, 9e-2, 'global sensitivity',color='black',fontsize=13,rotation=0, zorder=5)
    
    return fig, ax, lightest_neutrino, effective_masses