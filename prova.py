import re
string = "UL2018_v2"
year = string
if isinstance(string,str):
    matches = re.findall(r"(?:20)?[0-3]\d",string)
    if matches:
      matches.sort(key=lambda x: len(x),reverse=True)
      year = int(matches[0])
      if year<2000:
        year += 2000

print(year)
