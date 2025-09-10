import pandas as pd
import os

# Import the data from the csv file
data = pd.read_csv('MWS_wire_data.csv')

#conversion factor
conv_factor = 1000 * .3048 # one foot is .3048 meters

#create ohms per meter column 
data['ohms_per_meter'] = round(data['RESISTANCE_NOM'] / conv_factor, 6)
#create lbs per meter column
data['lbs_per_meter'] = round(data['POUNDS_PER_1000FT'] / conv_factor, 8)

# Save updated data to CSV (overwrite the original file)
data.to_csv('MWS_wire_data.csv', index=False)

print("Wire data updated successfully!")
print(f"Updated {len(data)} wire gauge entries")
print("Fixed lbs_per_meter precision for thin wires (AWG 43-52)")
#save data to a pandas dataframe .py file
