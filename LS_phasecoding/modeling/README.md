
Install docker: https://docs.docker.com/install/

once you have docker installed:

`docker pull davidtingley/papers`

`docker run -p 8888:8888 -i -t --net=host davidtingley/papers /bin/bash`

`jupyter notebook`

copy & paste the output into your browser, it should look something like this: 

`http://localhost:8888/?token=0e5fa29efdff8b412dd906e222d9ea1414c10224ce5df082&token=0e5fa29efdff8b412dd906e222d9ea1414c10224ce5df082`

From your jupyter notebook you should be able to open/run any of the *.ipynb notebooks used in the paper












