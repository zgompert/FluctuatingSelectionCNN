#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('pwd', '')


# In[2]:


get_ipython().run_line_magic('cd', '/uufs/chpc.utah.edu/common/home/gompert-group4/projects/fluctCNN/CNNTrainingData')


# In[45]:


import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow import keras
from keras.preprocessing.image import ImageDataGenerator
from keras.utils import np_utils
from tensorflow.keras import layers


# In[4]:


import os
directory_path ='/uufs/chpc.utah.edu/common/home/gompert-group4/projects/fluctCNN/CNNTrainingData'


# In[5]:


# making a list of all csv files
csv_files = [file for file in os.listdir(directory_path) if file.endswith('.csv')]


# In[6]:


print(len(csv_files))
#print(csv_files)
df=pd.read_csv(csv_files[1], header=None, sep=' ', dtype='float32')


# In[7]:


# extracting labels from file name of each CSV and storing it in a vector "vec"
# my_list = csv_files[:]
# vec = [item[:5] for item in my_list]


# In[8]:


# normalizing value of environment in each even columns in each CSV file by using "min-max scaling"
# for file in csv_files:
#     filepath = file  # You can modify the filepath if the CSV files are in a different directory
#     filename = os.path.splitext(file)[0]  # Extract the filename without extension

#     df = pd.read_csv(filepath)  # Read the CSV file into a pandas DataFrame

#     # Normalize even columns using min-max scaling
#     even_cols = df.columns[0::2]  # Select even columns
#     df[even_cols] = ((df[even_cols] + 2)/ 4)

#     # Save the normalized DataFrame to a new CSV file with the same name
#     normalized_filepath = f"{filename}_normalized.csv"
#     df.to_csv(normalized_filepath, index=False)

#     print(f"Normalized file saved: {normalized_filepath}")


# In[96]:


# csvn_files = [file for file in os.listdir(directory_path) if file.endswith('_normalized.csv')]
#print(len(csv_files))


# In[97]:


# a =[]
# for matrix in csv_files:
#     df=pd.read_csv(matrix, header=None, sep=' ', dtype='float32')
#     if df.shape==(10, 20):
#         a.append(1)
# print(sum(a))


# In[11]:


# Load the CSV matrices into numpy arrays
matrices = []
out_vec = []
for file in csv_files:
    df = pd.read_csv(file, header=None, sep=' ', dtype='float32')
    #df = df.iloc[1:, 1:].values
    #df = df.astype('float32')
    matrix = df
    matrices.append(matrix)
    # appending first five letters to out_vec as labels
    label = file[:5]
    out_vec.append(label)

#print(matrices)
# Convert the list of matrices into a single numpy array
# input_data = np.array(matrices)
# output_data = np.array(out_vec)


# In[65]:


# # to check how many entries are there in the matrices object
# a =[]
# for matrix in matrices:
#     a.append(1)
# print(sum(a))


# a =[]
# for matrix in out_vec:
#     a.append(1)
# print(sum(a))


# In[66]:


# # To double check whether the shapes of the matrices are same or not
# shapes=[matrix.shape for matrix in matrices]
# if len(set(shapes))>1:
#     print("matrices have different shapes")
# else:
#     print("matrices have same sahpe")


# In[ ]:


# converting input arrays of matrices into 
input_arrays = np.stack(matrices)
#print(input_arrays)
input_arrays.shape


# In[17]:


# converting label arrays of matrices into 
# import keras as keras
# from keras_utils import to_categorical
# out_vecc=to_categorical(out_vec)
# output_arrays = np.stack(out_vecc)
# #print(input_arrays)
# output_arrays.shape


# In[67]:


labels =out_vec
label_mapping = {label: index for index, label in enumerate(set(labels))}
print(label_mapping)


# In[21]:


integer_labels = np.array([label_mapping[label] for label in labels])


# In[69]:


# print(integer_labels)
categorical_labels = tf.one_hot(integer_labels, depth=len(label_mapping))


# In[70]:


out_labels = np.stack(categorical_labels)
out_labels.shape


# In[71]:





# In[25]:


# converting numpy arrays to tf tensors
input_tensor = tf.convert_to_tensor(input_arrays)
# converting output arrays to tf tensors
Label_tensor = tf.convert_to_tensor(out_labels)
print(tf.is_tensor(Label_tensor))
# to know shapes of input and label tensors
input_tensor.shape 
Label_tensor.shape


# In[36]:


input_arrays.shape
input_arr = np.reshape(input_arrays, (10000, 10, 20, 1))


# In[73]:


# # making values to two decimal places
# input_tensorr = tf.round(tf.multiply(input_tensor, 1000))/1000
# input_tensorr


# In[72]:


# # reshaping input tensor
# input_tensorrr = tf.reshape(input_tensorr, shape=(10000, 10, 20, 1))
# #input_tensorrr = tf.expand_dims(input_tensorr, axis=-1)
# input_tensorrr.shape
# #input_tensorrrx[[1]]


# In[74]:


## dividing input and label data into training, testingsets 4k, 1k size
from sklearn.model_selection import train_test_split
input_train, input_test, output_train, output_test = train_test_split(input_arr, out_labels, test_size=0.1, random_state=42)


# In[78]:


#print(output_test)


# In[94]:


model = keras.Sequential()
model.add(layers.Conv2D(32, kernel_size=(2, 2), activation="relu", input_shape=(10, 20, 1)))
model.add(layers.MaxPooling2D(pool_size=(1, 1)))
model.add(layers.Conv2D(64, kernel_size=(2, 2), activation="relu"))
model.add(layers.MaxPooling2D(pool_size=(2, 2)))
model.add(layers.Flatten())
model.add(layers.Dense(64, activation="relu"))
model.add(layers.Dense(2, activation="sigmoid"))

optimizer = keras.optimizers.Adam(learning_rate=0.001)
# Compile the model
model.compile(optimizer=optimizer, loss="categorical_crossentropy", metrics=["accuracy"])

# Train the model
model.fit(input_train, output_train, batch_size=64, epochs=25, validation_split=0.1)

# Evaluate the model on the test set
loss, accuracy = model.evaluate(input_test, output_test)
# print("Test Loss:", loss)
# print("Test Accuracy:", accuracy)

