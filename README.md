#Installation guide

Option 1: Docker

Requirements: 
-Docker client

Pull docker immage
```docker pull chrisamson/shiny-popoff-app:latest```

Run docker instance
```docker run chrisamson/shiny-popoff-app:latest```

Open port 3838 of the docker instance in browser to access PopOff

Option 2: Docker

Requirements:
-R
-plyr >= 1.8.9
-purr >= 1.0.1
-RSQLite >= 2.3.0
-scales >= 1.3.0
-shiny >= 1.8.0
-stringi >= 1.7.12
-tidyr >= 1.3.0
-DT >= 0.31
-ggplot2 >= 3.4.4
-ggpubr >= 0.6.0
-dplyr >= 1.1.0
-grid >= 4.2.2
-gridExtra >= 2.3

Clone github repository
```gh repo clone chrisamson/PopOff```

Change direcorty to 
```PopOff/POP_OFF```

Run app.R

Open http://127.0.0.1:6973/ in browser to acces PopOff
