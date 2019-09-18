Adding the folder 'youramrvacdir/tools/python' to the Python path:

Windows users: Start -> control panel -> 'edit the system environment variables' -> tab 'advanced' -> 'environment variables'
	       search for PYTHONPATH (or create one if it is not already there) and add 'youramrvacdir/tools/python'.

Mac & Linux users can add the following line to their PYTHONPATH in ~/.bashprofile (Mac) or ~/.bashrc (Linux)

                   export PYTHONPATH="${PYTHONPATH}:../.."

where you replace '../..' with the actual path to 'youramrvacdir/tools/python'.

--------------------------------------------------------------

Users that not want to modify their PYTHONPATH, can add these two lines at the beginning of every Python script where the tools are to be used:

                   import sys
                   sys.path.append('path to amrvac_tools')

where you replace 'path to amrvac_tools' with the actual path 'youramrvacdir/tools/python'.