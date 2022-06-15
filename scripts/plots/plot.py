
import matplotlib.pyplot as plt
import pandas as pd
  
  
# Initialize the lists for X and Y
data = pd.read_csv('/Users/ant/Documents/con/coocsr.csv')
  
df = pd.DataFrame(data)
print(df)

df['sparskit'] = df['sparskit'] /df['ours']
df['mkl'] = df['mkl'] /df['ours']
df['taco'] = df['taco'] /df['ours']
df['ours'] = df['ours'] /df['ours']

df.plot(x="matrix", y=['sparskit','mkl','taco'], kind="bar")


#df.groupby(['matrix']).plot.bar(x = 'matrix', y = ['sparskit','mkl','taco'])
  
# # Plot the data using bar() method
# plt.bar(X, Y, color='g')
# plt.title("Students over 11 Years")
# plt.xlabel("Years")
# plt.ylabel("Number of Students")
  
# # Show the plot
# plt.show()