from airflow import DAG
from airflow.operators.python import PythonOperator
from airflow.utils.dates import days_ago
from airflow.models import Variable
from airflow.providers.postgres.hooks.postgres import PostgresHook
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import boto3
from io import BytesIO
from datetime import datetime

S3_BUCKET = Variable.get("S3_BUCKET")
AWS_ACCESS_KEY_ID = Variable.get("AWS_ACCESS_KEY_ID")
AWS_SECRET_ACCESS_KEY = Variable.get("AWS_SECRET_ACCESS_KEY")
AWS_REGION = Variable.get("AWS_REGION") 


default_args = {
    'owner': 'airflow',
    'retries': 1,
    'retry_delay': 300,
}

dag = DAG(
    'smiles_to_s3',
    default_args=default_args,
    description='ETL for SMILES data, transform molecular properties, and upload to S3',
    schedule_interval='@daily',
    start_date=days_ago(1),
    catchup=False,
)


def extract_smiles_data():
    postgres_hook = PostgresHook(postgres_conn_id = 'postgres_smiles_db')
    engine = postgres_hook.get_sqlalchemy_engine()
    query = """
    SELECT smiles FROM molecules;
    """
    data = pd.read_sql(query, engine)
    return data.to_dict('records')


def transform_data(**kwargs):
    data = kwargs['ti'].xcom_pull(task_ids='extract_data')
    df = pd.DataFrame(data)
    
    
    def calc_properties(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return {
                'MW': Descriptors.MolWt(mol),
                'logP': Descriptors.MolLogP(mol),
                'TPSA': Descriptors.TPSA(mol),
                'H_Donors': Descriptors.NumHDonors(mol),
                'H_Acceptors': Descriptors.NumHAcceptors(mol),
                'Lipinski': (Descriptors.MolWt(mol) < 500 and Descriptors.MolLogP(mol) < 5)
            }
        return {}

    properties = df['smiles'].apply(calc_properties)
    properties_df = pd.DataFrame(properties.tolist())
    
    result = pd.concat([df, properties_df], axis=1)
    return result.to_dict('records')


def save_to_s3(**kwargs):
    data = kwargs['ti'].xcom_pull(task_ids='transform_data')
    df = pd.DataFrame(data)
    
    
    s3_client = boto3.client('s3',
                             aws_access_key_id=AWS_ACCESS_KEY_ID,
                             aws_secret_access_key=AWS_SECRET_ACCESS_KEY,
                             region_name=AWS_REGION)

    with BytesIO() as output:
        df.to_excel(output, index=False)
        output.seek(0)
        s3_client.put_object(Bucket=S3_BUCKET, Key=f"smiles_data_{datetime.now().date()}.xlsx", Body=output.getvalue())

extract_data = PythonOperator(
    task_id='extract_data',
    python_callable=extract_smiles_data,
    dag=dag,
)

transform_data = PythonOperator(
    task_id='transform_data',
    python_callable=transform_data,
    provide_context=True,
    dag=dag,
)

save_to_s3 = PythonOperator(
    task_id='save_to_s3',
    python_callable=save_to_s3,
    provide_context=True,
    dag=dag,
)

extract_data >> transform_data >> save_to_s3
