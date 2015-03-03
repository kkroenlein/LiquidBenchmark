import pandas as pd
import bs4
import re

reg = re.compile(r"""<td align="center">\n""")

def parse_tr(tr):
    bits = re.findall("\>(.*?)\<", str(tr))
    data = []
    print(bits)
    for bit in bits:
        if "\xc2\xb11" in bit:
            bit = bit.replace("\xc2\xb11", "")
        try:
            value = float(bit)
            data.append(value)
        except Exception as e:
            print(e)
    print(data)
    return data
    

def parse_page(soup):
    tables = soup.findAll("table")
    data = {}
    for table in tables:
        for key in ["Static dielectric", "Liquid density"]:
            if key not in data:
                if key in table.prettify() and "Physical properties" not in table.prettify():
                    trs = table.findAll("tr")
                    for tr in trs:
                        if key in tr.prettify():
                            print(tr)
                            data[key] = parse_tr(tr)
    return data


metadata = pd.read_table("./virtualchemistry.tab")
data = []

for cas in metadata.cas:
    print(cas)
    filename = "./pages/page_%s.html" % cas
    f = open(filename).read()
    soup = bs4.BeautifulSoup(f)
    current_data = parse_page(soup)
    print(current_data)
    if "Liquid density" in current_data:
        print(current_data["Liquid density"])
        try:
            density = current_data["Liquid density"][1]    
        except IndexError:
            density = np.nan
    if "Static dielectric" in current_data and len(current_data["Static dielectric"]) == 4:
        temperature, expt, gaff, opls = current_data["Static dielectric"]
        data.append(dict(cas=cas, density=density, temperature=temperature, expt=expt, gaff=gaff, opls=opls))

        

data = pd.DataFrame(data)
data
data.to_csv("./vchem.csv")
