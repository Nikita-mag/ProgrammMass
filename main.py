import requests
import json
import scipy.stats
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)

def check_spectrum(df):
    if len(df) == 0:
        return np.zeros(0, dtype=np.float64).reshape(-1, 2)
    df = np.asarray(df, dtype=np.float64, order="C")
    sz = df.shape
    kl = sz[1]
    if kl != 2:
                raise RuntimeError("Ошибка во входном формате спектра!")
    df = pd.DataFrame(columns=['m/z', 'I'], data=df)
    return df

def save_from_www(link):
    filename = link.split('/')[-1]
    print(filename)
    r = requests.get(link, allow_redirects=True)
    open(filename,"wb").write(r.content)

link1 = 'https://raw.githubusercontent.com/Nikita-mag/ProgrammMass/main/1.csv'
link2 = 'https://raw.githubusercontent.com/Nikita-mag/ProgrammMass/main/2.csv'
link3 = 'https://raw.githubusercontent.com/Nikita-mag/ProgrammMass/main/3.csv'

save_from_www(link1)
save_from_www(link2)
save_from_www(link3)


df =pd.read_csv('2.csv', sep=';')
df = check_spectrum(df)

#Creating a result dataframe
df_m = pd.DataFrame(columns=['m/z', 'Spectral entropy', 'RT', 'Identify', 'Height'])

#m/z
str1 = df['m/z'][df['I'] == df['I'].max()]
df_m = pd.DataFrame(str1)

#Spectral entropy
print('-'* 30)
spectral_entropy = scipy.stats.entropy(df)
spectral_entropy = spectral_entropy[spectral_entropy.shape[0]-1]
print("Spectral entropy is {}.".format(spectral_entropy))
df_m.insert (loc= len(df_m.columns) , column='Spectral entropy', value=format(spectral_entropy))

#Identify
spectral_entropy = str(spectral_entropy)
with open ('MoNA-export-LC-MS_Spectra.json', 'r') as file:
    data = json.load(file)
    for i in range (0,2991):
        a = data[i].get('compound', 0)
        b = data[i].get('metaData', 0)  # значение спектральной энтропии 15
        b1 = b[15].get('value', 0)
        b2 = b[9].get('value', 0)
        a1 = a[0].get('names', 0)
        a2 = a1[0].get('name', 0)
        if spectral_entropy == b1:
            df_m.insert(loc=len(df_m.columns), column='RT', value=b2)
            df_m.insert(loc=len(df_m.columns), column='Identify', value=a2)
            break
    else:
        print('Не найдена спектральная энтропия!')


#Height
str2 = df['I'].max()
df_m.insert (loc= len(df_m.columns) , column='Height', value=format(str2, '.2e'))

#Формирование и вывод спектрограммы
x = df["m/z"]
y = df["I"]
fig, ax = plt.subplots(figsize=(8, 8))
margins = {
    "left"   : 0.140,
    "bottom" : 0.150,
    "right"  : 0.990,
    "top"    : 0.920
}
fig.subplots_adjust(**margins)
ax.set_title("Result", fontsize=16)
ax.set_xlabel("m/z", fontsize=14)
ax.set_ylabel("Intensity", fontsize=14)
ax.stem(x, y, linefmt="-", markerfmt="none",orientation='vertical')

ax.grid(which="major", linestyle="--", color="gray", linewidth=0.7)
ax.grid(which="minor", linestyle="--", color="gray", linewidth=0.7)

ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.xaxis.set_minor_formatter(FormatStrFormatter("%.4f"))
ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.e'))

ax.tick_params(which='major', length=10, width=2)
ax.tick_params(which='minor', length=5, width=1, rotation = 90)
ax.grid()
plt.show()


