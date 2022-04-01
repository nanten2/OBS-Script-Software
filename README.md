# OBS-Script-Software

Software to create OBS files for the NANTEN2 telescope

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
