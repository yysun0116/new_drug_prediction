FROM mcs07/rdkit:2020.03.2

# Install ubuntu packages
RUN apt update \
 && apt-get install -y git \
 && apt-get install -y unzip \
 && apt-get install -y  wget \
 && apt-get install -y  curl \
 && apt-get install -y  time 

# Install python packages
RUN pip3 install --upgrade pip
RUN pip3 install cmake --upgrade
RUN pip3 install pandas==1.5.3 numpy==1.24.1 scikit-learn==1.2.0

# Download dependent files for scDrug
RUN git clone https://github.com/yysun0116/new_drug_prediction.git /scDrug


# Load shell script
COPY run.sh /app/
WORKDIR /app

CMD ["sh", "run.sh"]
