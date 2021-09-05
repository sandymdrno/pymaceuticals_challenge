#!/usr/bin/env python
# coding: utf-8

# ## Observations and Insights 

# 

# In[124]:


# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st


# In[103]:


# Study data files
mouse_metadata_path = "data/Mouse_metadata.csv"
study_results_path = "data/Study_results.csv"

# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)


# In[104]:


# Combine the data into a single dataset
# Display the data table for preview
combined_df = pd.merge(mouse_metadata, study_results, how="inner", on="Mouse ID")
combined_df


# In[105]:


# Checking the number of mice.
mouse_count = combined_df["Mouse ID"].count()
mouse_count


# In[106]:


# Getting the duplicate mice by ID number that shows up for Mouse ID and Timepoint. 
duplicate_df = combined_df[combined_df.duplicated(['Mouse ID', 'Timepoint'])]
duplicate_df


# In[107]:


# Optional: Get all the data for the duplicate mouse ID. 


# In[108]:


# Create a clean DataFrame by dropping the duplicate mouse by its ID.
clean_df = combined_df.drop_duplicates("Mouse ID")

clean_df


# In[109]:


# Checking the number of mice in the clean DataFrame.
mouse_count_new = clean_df["Mouse ID"].count()
mouse_count_new

#clean_df["Mouse ID"]


# ## Summary Statistics

# In[110]:


mice_df = pd.merge(clean_df,study_results,how='outer', on = "Mouse ID")
mice_df

del mice_df['Timepoint_x']
del mice_df['Metastatic Sites_x']
del mice_df['Tumor Volume (mm3)_x']            
mice_df

mice_df = mice_df.rename(columns={'Timepoint_y': 'Timepoint', 'Tumor Volume (mm3)_y': 'Tumor Volume',
                        'Metastatic Sites_y': 'Metastatic Site'})
mice_df


# In[111]:


# Groupby Mice DataFrame and Drug Regimen
mice_groupby_df = mice_df.groupby(["Drug Regimen"])
mice_groupby_df.head()


# In[112]:


# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen
tumor_mean = round(mice_groupby_df['Tumor Volume'].mean(),2)
tumor_median = round(mice_groupby_df['Tumor Volume'].median(),2)
tumor_variance = round(mice_groupby_df['Tumor Volume'].var(),2)
tumor_standarddev = round(mice_groupby_df['Tumor Volume'].std(),2)
tumor_sem = round(mice_groupby_df['Tumor Volume'].sem(),2)


summary_tumor_volume_df = pd.DataFrame({'Average Tumor Volume': tumor_mean,
                                        'Median Tumor Volume': tumor_median,
                                        'Variance Tumor Volume': tumor_variance,
                                        'Std of Tumor Volume': tumor_standarddev,
                                        'SEM of Tumor Volume': tumor_sem})

summary_tumor_volume_df
# Use groupby and summary statistical methods to calculate the following properties of each drug regimen: 
# mean, median, variance, standard deviation, and SEM of the tumor volume. 
# Assemble the resulting series into a single summary dataframe.


# In[113]:


# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen
# Using the aggregation method, produce the same summary statistics in a single line

summary_tbl = mice_groupby_df.agg(['mean','median','var','std','sem'])['Tumor Volume']
summary_tbl


# ## Bar and Pie Charts

# In[114]:


# Generate a bar plot showing the total number of timepoints for all mice tested for each drug regimen using Pandas.

total_timepoints = mice_groupby_df.count().reset_index()
drugs = total_timepoints[["Drug Regimen","Mouse ID"]]
drugs = drugs.set_index("Drug Regimen")
drugs.plot(kind="bar", legend = False, rot = 50,ylabel = '# of mice')

plt.title("Number of mice tested with each drug")
plt.show

plt.show()


# In[115]:


# Generate a bar plot showing the total number of timepoints for all mice tested for each drug regimen using pyplot.

total_timepoints = mice_groupby_df.count()['Timepoint']
plt.ylabel('Number of Timepoints')
plt.title('Timepoints of all mice tested')
total_timepoints.plot.bar(legend=False,rot=50)
plt.show()


# In[116]:


# Generate a pie plot showing the distribution of female versus male mice using Pandas
# Collect vaule count of sex

sex_count = mice_df["Sex"].value_counts()
sex_count

sex_count.plot(kind="pie", autopct="%1.1f%%")
plt.show()


# In[184]:


# sex_count = mice_df.groupby(["Mouse ID", "Sex"]).count().reset_index()[["Mouse ID", "Sex"]]["Sex"].value_counts()

sex_count = mice_df.groupby(["Mouse ID", "Sex"]).count()
sex_count = sex_count.reset_index()
sex_count = sex_count[["Mouse ID", "Sex"]]

mouse_sex_df = sex_count["Sex"].value_counts()

mouse_sex_df = mouse_sex_df.to_frame()
mouse_sex_df


# In[128]:


# # Generate a pie plot showing the distribution of female versus male mice using Pyplot
# #mice_sex = pd.DataFrame
# sex_count = mice_df["Sex"].value_counts()
# mouse_sex_df = pd.DataFrame(sex_count)
# mouse_sex_df


# In[181]:


get_ipython().run_line_magic('matplotlib', 'notebook')


# In[182]:


# Generate a pie plot showing the distribution of female versus male mice using pyplot
# mice_sex = mice_groupby_df.plot(kind="pie", y=sex_list, title=("Sex of Mice " + mouse_id))
mice_sex = mouse_sex_df.plot(kind="pie", subplots=True, title=("Sex of Mice"), autopct = "%1.1f%%")

plt.show()
plt.axis("equal")
plt.tight_layout()


# ## Quartiles, Outliers and Boxplots

# In[145]:


# Calculate the final tumor volume of each mouse across the treatment regimen Capomulin
capomulin_vol = combined_df.loc[combined_df["Drug Regimen"] == "Capomulin"]

# Start by getting the last (greatest) timepoint for each mouse
capomulin_last = capomulin_vol.groupby("Mouse ID").max()["Timepoint"]
capomulin_vol = pd.DataFrame(capomulin_last)

capomulin_last
capomulin_vol


# In[146]:


# Calculate the final tumor volume of each mouse across the treatment regimen Ramicane
ramicane_vol = combined_df.loc[combined_df["Drug Regimen"] == "Ramicane"]

# Start by getting the last (greatest) timepoint for each mouse
ramicane_last = ramicane_vol.groupby("Mouse ID").max()["Timepoint"]
ramicane_vol = pd.DataFrame(ramicane_last)

ramicane_last
ramicane_vol


# In[148]:


# Calculate the final tumor volume of each mouse across the treatment regimen Infubinol
infubinol_vol = combined_df.loc[combined_df["Drug Regimen"] == "Infubinol"]

# Start by getting the last (greatest) timepoint for each mouse
infubinol_last = infubinol_vol.groupby("Mouse ID").max()["Timepoint"]
infubinol_vol = pd.DataFrame(infubinol_last)

infubinol_last
infubinol_vol


# In[144]:


# Calculate the final tumor volume of each mouse across the treatment regimen Ceftamin
ceftamin_vol = combined_df.loc[combined_df["Drug Regimen"] == "Ceftamin"]

# Start by getting the last (greatest) timepoint for each mouse
ceftamin_last = ceftamin_vol.groupby("Mouse ID").max()["Timepoint"] 
ceftamin_vol = pd.DataFrame(ceftamin_last)
ceftamin_last
ceftamin_vol 


# In[193]:


capomulin_vol.columns = ["capomulin"]
ramicane_vol.columns = ["ramicane"]
ceftamin_vol.columns = ["ceftamin"]
infubinol_vol.columns = ["infubinol"]
dfs = [
    capomulin_vol,
    ramicane_vol,
    ceftamin_vol,
    infubinol_vol
]
all_drugs = pd.concat(dfs, axis=1)
display(all_drugs.head(50))
display(all_drugs.tail(50))


# In[168]:


# Merge this group df with the original dataframe to get the tumor volume at the last timepoint
grp1_drugs = pd.merge(capomulin_vol,ramicane_vol, on = "Mouse ID", how = "outer",suffixes=('_capomulin', '_ramicane'))
tumor_vol_last_timepoint1 =  pd.merge(clean_df, grp1_drugs, on = "Mouse ID", how = "outer",suffixes=('_clean_df', '_grp1'))
grp2_drugs = pd.merge(ceftamin_vol,infubinol_vol, on = "Mouse ID", how = "outer",suffixes=('_ceftamin', '_infubinol'))
tumor_vol_last_timepoint =  pd.merge(tumor_vol_last_timepoint1, grp2_drugs, on = "Mouse ID", how = "outer",suffixes=('_clean', '_grp2'))
# drug_combined = pd.merge(grp1_drugs, grp2_drugs, on ="Mouse ID", how = "outer", suffixes=('_capomulin', '_ramicane','_ceftamin', '_infubinol'))
drug_combined = pd.merge(grp1_drugs,grp2_drugs, on ="Mouse ID", how = "outer", suffixes=('_capomulin', '_ramicane','_ceftamin', '_infubinol'))
                                                                                    
display (tumor_vol_last_timepoint)
display (drug_combined)                                                                                   


# In[194]:


# Put treatments into a list for for loop (and later for plot labels)
treatment_list = ["Capomulin", "Ramicane","Ceftamin", "Infubinol"]
treatment_list


# In[ ]:


# Create empty list to fill with tumor vol data (for plotting)
tumor_vol_list = []

# Calculate the IQR and quantitatively determine if there are any potential outliers. 
drug_combined
    
    # Locate the rows which contain mice on each drug and get the tumor volumes
    
    
    # add subset 
    
    
    # Determine outliers using upper and lower bounds
    


# In[ ]:


# Generate a box plot of the final tumor volume of each mouse across four regimens of interest


# ## Line and Scatter Plots

# In[ ]:


# Generate a line plot of tumor volume vs. time point for a mouse treated with Capomulin


# In[ ]:


# Generate a scatter plot of average tumor volume vs. mouse weight for the Capomulin regimen


# ## Correlation and Regression

# In[ ]:


# Calculate the correlation coefficient and linear regression model 
# for mouse weight and average tumor volume for the Capomulin regimen


# In[ ]:




