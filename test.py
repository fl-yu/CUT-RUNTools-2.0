



import json 




with open('/Users/fulongyu/github/CUT-RUNTools-2.0/sc-config.json') as f:
  data = json.load(f)

# Output: {'name': 'Bob', 'languages': ['English', 'Fench']}
print(data)


