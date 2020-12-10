### Load and Run CNCDriver Container Instructions

#### 1. Download the docker desktop:
- Visit the docker official site at https://www.docker.com/get-started
- Download the docker desktop for your operating system


#### 2. Load the docker image:

- Extract the image:

    ```bash
    tar xvz cnc_driver_docker.tar.gz
    ``` 
    
- Load the docker image from the tar archive using the following command:

    - For macOS and Linux:
    
        ```bash
        docker load < cncdriver.tar
        ```   
    - For Windows:
        ```bash
        docker load -i cncdriver.tar
        ``` 
- NOTE: Depending on your system privilege level, you may need to prefix this command with ```sudo```


#### 4. Modify the paths in the shell script:
-  Open ```run_docker_cncdriver.sh``` script and modify the paths to your local paths where they indicated so before you attempt to run the script

-  You will also need to specifiy the container's paths in your own R script. A sample script ```sample_cncdriver.R```  was provided for you to serve as tempate or example

- Finally, run the following command:

    ```bash
    ./run_docker_cncdriver.sh
    ```
