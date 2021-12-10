from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image
import moviepy.editor as mp

names = ['A', 'CA', 'CH', 'D', 'F', 'Ma', 'Mr'] 
labels = ['Antigen concentration', 'Anti-inflammatory Cytokine',
          'pro-inflammatory Cytokine', 'Tissue Damage', 
          'Antibodies', 'Activated Antigen Presenting Cells',
          'Resting Antigen Presenting Cells']

for i in range(7):

    # Generate each frame
    for n in range(0,69888,416):
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        data = pd.read_csv(f"output/{names[i]}_{n}.csv", header=None, names =('x','y','z','value'))
        #print(data.head())
        data = data.loc[(data['y']) >= 5]
        ax.scatter(data.x, data.y, data.z, c=data.value, cmap='viridis', 
                label=f"{labels[i]}", linewidth=0.5)
        #ax.set_ylim3d(5, 10)
        fig.colorbar(plt.cm.ScalarMappable(cmap = 'viridis'), ax = ax)
        ax.legend()
        plt.savefig(f"output/{names[i]}_{n}.png")
        plt.close()
        
    images = [Image.open(f"output/{names[i]}_{n}.png") for n in range(0,69888,416)]

    images[0].save(f"output/animations/{names[i]}.gif", save_all=True, append_images=images[1:], duration=100, loop=0)
    
    clip = mp.VideoFileClip(f"output/{names[i]}.gif")
    clip.write_videofile(f"output/{names[i]}.mp4")