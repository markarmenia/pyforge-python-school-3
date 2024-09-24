def lambda_handler(event, context):
    names = event.get('names', ['World'])  # This line is correct
    greetings = [f'Hello, {name}!' for name in names]
    return {
        'statusCode': 200,
        'body': greetings
    }

import boto3
import json

client = boto3.client('lambda', region_name='us-east-1')

event = {
    "names": ["Alice", "Bob", "Charlie"]  # Changed 'name' to 'names'
}

response = client.invoke(
    FunctionName='HelloStudentFunction',  # Replace with your Lambda function's name
    InvocationType='RequestResponse',
    Payload=json.dumps(event),
)

response_payload = json.loads(response['Payload'].read())
print(response_payload)
