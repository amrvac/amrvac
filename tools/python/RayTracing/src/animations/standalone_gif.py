"""
Standalone script to generate a .gif file from a pack of user-specified .png files.
Can run separately, change containing folders accordingly.
"""

from animations import animate_series

halpha_name  = "k6pi_05MK_Halpha_x"
faraday_name = "k6pi_05MK_Faraday_x"

if __name__ == '__main__':

    filenames_halpha  = []
    filenames_faraday = []
    for i in range(93):
        filenames_halpha.append("png_files/halpha/{}_{}.png".format(halpha_name, str(i).zfill(4)))
        filenames_faraday.append("animations/png_files/faraday/{}_{}.png".format(faraday_name, str(i).zfill(4)))

    animate_series.create_gif(filenames_halpha)
    animate_series.create_gif(filenames_faraday)

