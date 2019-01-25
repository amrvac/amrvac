"""
Module to animate a series of .dat files.
Depending on the user choice, a movie can be made showing an evolution in H-alpha or others.
The appropriate conditions must be specified in settings.py: the path to the folder where the .dat files are,
together with the starting and stopping index of the snapshots.

Created on 22 Jan 2019

@author: Niels Claes
"""
import settings
import os, sys, glob
import matplotlib.pyplot as plt
from views import h_alpha_view, faraday_view
from views import plotting
from dataIO.reduce_data import ProcessData
import print_tools
import imageio


def check_settings():
    """
    Checks folders and settings before program execution.
    """

    print("Saving animations...")
    if not os.path.isdir(settings.folderpath):
        sys.exit("Folder not found: {}".format(settings.folderpath))
    if settings.animate_end < settings.animate_start:
        sys.exit("Start animation ({}) has larger index than end ({})".format(settings.animate_start, settings.animate_end))
    if not os.path.isdir("animations/png_files"):
        os.mkdir("animations/png_files")
        print("Created directory animations/png_files")
    if not os.path.isdir("animations/gif_files"):
        os.mkdir("animations/gif_files")
        print("Created directory animations/gif_files")
    if not os.path.isdir("animations/png_files/halpha"):
        os.mkdir("animations/png_files/halpha")
        print("Created directory animations/png_files/halpha")
    if not os.path.isdir("animations/png_files/faraday"):
        os.mkdir("animations/png_files/faraday")
        print("Created directory animations/png_files/faraday")
        print("\n")
    if not (settings.line_of_sight == "x" or settings.line_of_sight == "y" or settings.line_of_sight == "z"):
        sys.exit("Line of sight is not specified correctly in settings. Allowed values are 'x', 'y', 'z'.\n"
                 "'all' is currently not supported (yet).\n"
                 "Current value: '{}'".format(settings.line_of_sight))
    return


def search_files():
    """
    Search for .dat files in the folder specified in settings.py
    :return: (List) Sorted list containing the .dat filenames as strings.
    :raise: IOError if no files are found.
    """
    main_dir = os.getcwd()
    os.chdir(settings.folderpath)
    files = []
    for filename in sorted(glob.glob('*.dat')):
        files.append(filename)

    if not files:
        raise IOError("No files found, folder empty?")

    os.chdir(main_dir)
    return files


def get_animation_sequence():
    """
    Returns the appropriate animation sequence as specified by start and end
    in settings.py.
    :return: (List) Slice of complete .dat file list, corresponding to required interval.
    """
    files = search_files()
    return files[settings.animate_start:settings.animate_end+1]


def check_existence(filename):
    """
    Checks if .png file already exists. If so and if overwriting is
    disabled, the current .dat file will be skipped.
    :param filename: (String) The filename.
    :return: True if .png file exists and if overwriting is disabled. False otherwise.
    """
    if settings.overwrite_png:
        return False
    if os.path.isfile(filename):
        return True
    return False


def create_gif(images, duration=None):
    """
    Creates a .gif file from a given list of image paths.
    :param images: Sorted list containig the filepaths to the desired images.
    :param duration: Duration of each frame, in seconds.
    """
    if duration is None:
        duration = settings.frame_duration

    imageArray = []
    print("\nGENERATING GIF FILE...")
    counter = 0
    for im in images:
        imageArray.append(imageio.imread(im))
        counter += 1
    print("Iterated over {} images".format(counter))
    gifname = "animations/gif_files/{}.gif".format(print_tools.trim_filename(images[0])[:-5])
    imageio.mimsave(gifname, imageArray, duration=duration)
    print(".gif file saved to {}".format(gifname))
    return


def check_to_remove_png(pngfiles):
    """
    Checks and asks again if png files must be deleted. This can not be undone.
    :param pngfiles: (List) List containing the paths to the .png files.
    """
    if settings.remove_pngfiles:
        print("Files >>{}.png<< through >>{}.png<< will be removed.".format(print_tools.trim_filename(pngfiles[0]),
                                                                            print_tools.trim_filename(pngfiles[-1])))
        usr_input = input("Continue? Y/N: ")
        if usr_input == "yes" or usr_input == "Y" or usr_input == "y":
            for png_file in pngfiles:
                os.remove(png_file)
            print("Files removed.\n")
        else:
            print("Files not removed.\n")
    return


def check_to_remove_interpolated():
    """
    Checks and asks again if interpolated .npy files must be deleted. This can not be undone.
    Assumes that files have default naming convention defined in physics/ionization.py, eg. ion_ and f_
    """
    anim_seq = get_animation_sequence()
    f_npy_files   = []
    ion_npy_files = []
    for filename in anim_seq:
        ion_npy_files.append("interpolated_files/ion_{}.npy".format(filename[:-4]))
        f_npy_files.append("interpolated_files/f_param_{}.npy".format(filename[:-4]))


    if settings.remove_npyfiles:
        print("Files >>{}.npy<< through >>{}.npy<< will be removed.".format(print_tools.trim_filename(ion_npy_files[0]),
                                                                            print_tools.trim_filename(ion_npy_files[-1])))
        print("Files >>{}.npy<< through >>{}.npy<< will be removed.".format(print_tools.trim_filename(f_npy_files[0]),
                                                                            print_tools.trim_filename(f_npy_files[-1])))
        usr_input = input("Continue? Y/N: ")
        if usr_input == "yes" or usr_input == "Y" or usr_input == "y":
            for npy_file in ion_npy_files:
                os.remove(npy_file)
            for npy_file in f_npy_files:
                os.remove(npy_file)
            print("Files removed.\n")
        else:
            print("Files not removed.\n")
    return



def generate_animation():
    """
    Main body of this module. Generates the animation and checks the appropriate things.
    """
    check_settings()
    anim_seq = get_animation_sequence()

    # get current source directory
    main_dir = os.getcwd()

    figlist_halpha = []
    figlist_faraday = []

    for filename in anim_seq:
        # Set new file name
        settings.filename = settings.folderpath + "/" + filename
        # Get file number
        file_nb   = filename[-8:-4]
        # Remove .dat extension and file number
        file_trim = filename[:-8]

        # Create save names for figures
        savename_halpha  = file_trim + "_Halpha_{}_{}".format(settings.line_of_sight, file_nb)
        figname_halpha   = "animations/png_files/halpha/{}.png".format(savename_halpha)
        savename_faraday = file_trim + "_Faraday_{}_{}".format(settings.line_of_sight, file_nb)
        figname_faraday  = "animations/png_files/faraday/{}.png".format(savename_faraday)
        # Append to arrays for .gif creation
        figlist_halpha.append("{}".format(figname_halpha))
        figlist_faraday.append("{}".format(figname_faraday))

        # Check if files already exist, prevent loading of .dat file if so
        if check_existence(figname_halpha) or not settings.halpha:
            halpha_skip = True
        else:
            halpha_skip = False
        if check_existence(figname_faraday) or not settings.faraday:
            faraday_skip = True
        else:
            faraday_skip = False

        if halpha_skip and faraday_skip:
            print("Both .png files already exist. Skipping {}".format(filename))
            continue


        # Read .dat files
        os.chdir(settings.folderpath)
        dat_file = open(filename, "rb")
        # Change back to source directory
        os.chdir(main_dir)
        data = ProcessData(dat_file)

        if settings.halpha:
            if halpha_skip:
                print(".png file already exists: {}".format(figname_halpha))
                pass

            cbtitle = "H-alpha intensity [cgs]"

            intensity = h_alpha_view.get_intensity(data, axis=settings.line_of_sight)

            title = "H-alpha view along {}".format(settings.line_of_sight) \
                    + "\n" + file_trim + "\n" + " {}".format(file_nb)
            plotting.plot_surface(data, intensity, title=title, line_of_sight=settings.line_of_sight,
                                  cmap=settings.cmap_halpha, logscale=settings.logscale_halpha, cblabel=cbtitle,
                                  cbformat=settings.cbar_format_halpha)
            plt.savefig(figname_halpha)
            plt.close()
            print("\n>>>Figure saved to {}".format(figname_halpha))


        if settings.faraday:
            if faraday_skip:
                print(".png file already exists: {}".format(figname_faraday))
                pass
            if settings.include_lambda:
                cbtitle = r"Faraday Rotation $\beta$ [rad]" + "\n$\lambda$ = {} cm".format(settings.faraday_lambda)
            else:
                cbtitle = r"Faraday RM [rad m$^{-2}$]"

            faraday = faraday_view.get_faraday_strength(data, axis=settings.line_of_sight)
            title = "Faraday view along {}".format(settings.line_of_sight) \
                    + "\n" + file_trim + "\n" + "{}".format(file_nb)
            plotting.plot_surface(data, faraday, title=title, line_of_sight=settings.line_of_sight,
                                  cmap=settings.cmap_faraday, logscale=settings.logscale_faraday,
                                  cblabel=cbtitle, cbformat=settings.cbar_format_faraday)
            plt.savefig(figname_faraday)
            plt.close()
            print(">>>Figure saved to {}".format(figname_faraday))

        print("\n======================================\n")

    if settings.halpha:
        create_gif(figlist_halpha)
    if settings.faraday:
        create_gif(figlist_faraday)

    # Check for removal
    if settings.halpha:
        check_to_remove_png(figlist_halpha)
        check_to_remove_interpolated()
    if settings.faraday:
        check_to_remove_png(figlist_faraday)
        check_to_remove_interpolated()

    return