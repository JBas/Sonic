# import cv2
# import re
# import numpy as np
# import glob

# out = cv2.VideoWriter("video_001.avi", cv2.VideoWriter_fourcc(*"DIVX"), 60, (512, 512))
# # print("Created out!")

# # img_array = []
# # files = glob.glob('frames/*.png')
# # files.sort(key=lambda x:[int(c) if c.isdigit() else c for c in re.split(r'(\d+)', x)])
# # for filename in files:
# #    img = cv2.imread(filename)
# #    img_array.append(img)

# # for i in range(len(img_array)):
# #    out.write(img_array[i])
# out.release()


import os, re, glob
import moviepy.video.io.ImageSequenceClip

fps = 60

files = glob.glob("frames/*.png")
files.sort(key=lambda x:[int(c) if c.isdigit() else c for c in re.split(r"(\d+)", x)])

clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(files, fps=fps)
clip.write_videofile("my_video.mp4")