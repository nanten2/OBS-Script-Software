# OBS-Script-Software

Software to create OBS files for the NANTEN2 telescope

## Installation

1. Download the file from release assets
2. **\[Ubuntu\]** Extract (unzip) the file  
   **\[MacOS\]** Double click the disk image, move it to Application folder
3. To launch the app,  
   **\[Ubuntu\]** hit a command `./path/to/extracted/obs_pS` (working hard to enable double clicking launching)  
   **\[MacOS\]** Click the obs_pS.app in Launchpad

- `obs2_ubuntu-18.04.zip` supports Ubuntu 18.04+
- `obs2_macos-latest.dmg` supports MacOS 10.15(Catalina)+

### Command Line Installation

- Ubuntu:

    ```shell
    # download and install the app:
    $ curl -L -O https://github.com/nanten2/OBS-Script-Software/releases/latest/download/obs2_ubuntu-18.04.zip
    $ unzip -q path/to/downloaded/obs2_<OperatingSystem>.zip
    # then launch it:
    $ ./path/to/unzipped/obs2_ubuntu-18.04/obs_pS  # on Ubuntu
    ```

- MacOS

    In preparation...

## For Developers

- Implementation is documented [here](https://obs-script-software.readthedocs.io/en/latest/).
- To debug, run the entry point script from console:
  
  ```shell
  poetry run python scripts/obs_pS.py
  ```

  or execute the application directly from console:

  ```shell
  ./obs_pS.app/Contents/MacOS/obs_pS
  ```

- To release updated version, use `git tag` command like below.

    ```shell
    # commit all changes for the new release
    $ git tag v1.0.0  # create a tag (the name must be v*.*.*)
    $ git push origin v1.0.0  # push the tag to the remote
    ```

    Then the application is automatically built and compiled, and uploaded [here](https://github.com/nanten2/OBS-Script-Software/releases/latest).
