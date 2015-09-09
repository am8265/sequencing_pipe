### Uses curl. Please make sure you have the module

# Set secret token specific to your REDCap project
TOKEN="3B7C39D9720D72EB8B03CD468B060C1A"

# Set the url to the api (ex. https://YOUR_REDCAP_INSTALLATION/api/)
SERVICE="https://wchredcap.cumc.columbia.edu/redcap/api"

# UPLOAD a file (/path/to/myfile.txt) to field MYFILE for the record MYRECORDNUM and 
# use 'display.txt' as the file name REDCap will display.
# Note the use of '@' to get curl to read the file as a file attachment

curl    --form token=${TOKEN} \
        --form content=file \
	--form action=import \
	--form field=MYFILE \
	--form record=MYRECORDNUM 
	--form file="~/display.txt" \
	${SERVICE} 
