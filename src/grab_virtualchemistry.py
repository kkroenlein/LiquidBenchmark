import time
import urllib
import pandas as pd

data = pd.read_table("./virtualchemistry.tab")


for cas in data.cas:
    print(cas)
    path = "http://virtualchemistry.org/molecules/%s/index.php" % cas
    filename = "./pages/page_%s.html" % cas
    urllib.urlretrieve(path, filename)
    time.sleep(0.33)
