import pandas as pd
from pymongo import MongoClient

# Load the data from the file
data = pd.read_csv('25.H_sapiens.goa', sep='\t',comment='!', header=None)
colnames = ['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'Qualifier', 'GO_ID', 'DB:Reference', 'Evidence_Code', 'With_or_From', 'Aspect', 'DB_Object_Name', 'DB_Object_Synonym', 'DB_Object_Type', 'Taxon', 'Date', 'Assigned_By', 'Annotation_Extension', 'Gene_Product_Form_ID']
data.columns = colnames
data

# Filter the data to get the interesting ones
data2 = data[data['Qualifier'] == 'involved_in']

# Group the data and delete the duplicates (if not this can cause problems in the similarity counts)
final_data = (
    data2.groupby('DB_Object_Symbol')['GO_ID']
    .apply(lambda x: list(set(x)))
    .reset_index()
)

# Rename the columns to make more intuitive
final_data.columns = ['gene_name', 'go_id']
final_data

# Connect to MongoDB and insert the data
client = MongoClient('mongodb://localhost:27017/')
db = client['saRNAs']
collection = db['go_annotations']

# Insert
records = final_data.to_dict(orient='records')
collection.insert_many(records)

print(f"Inserted {len(records)} records into MongoDB.") 
client.close()