from plotly import graph_objects as go
import pandas as pd

# Load your updated data
data = pd.read_excel("20mer_hydrolysis_nested_algorithm_output.xlsx")  # Replace with your file path

# Create the 3D scatter plot
fig = go.Figure()

# Add scatter plot with lines and markers
fig.add_trace(go.Scatter3d(
    x=data['apex_rt'],  # Retention Time (min)
    y=data['monoisotopic_mass'],  # Mass (Da)
    z=data['sum_intensity'],  # Intensity
    mode='lines+markers',
    marker=dict(size=4, color='blue'),  # Marker size and color
    line=dict(width=2, color='blue')  # Line width and color
))

# Customize the layout with grey grid and white background
fig.update_layout(
    scene=dict(
        xaxis=dict(
            title='Retention Time (min)',
            showgrid=True,
            gridcolor='grey',  # Set gridlines color to grey
            zerolinecolor='grey'  # Ensure zero-line matches the gridlines
        ),
        yaxis=dict(
            title='Mass (Da)',
            showgrid=True,
            gridcolor='grey',  # Set gridlines color to grey
            zerolinecolor='grey'  # Ensure zero-line matches the gridlines
        ),
        zaxis=dict(
            title='Intensity',
            showgrid=True,
            gridcolor='grey',  # Set gridlines color to grey
            zerolinecolor='grey'  # Ensure zero-line matches the gridlines
        ),
        bgcolor='white'  # Overall scene background
    ),
    showlegend=True  # Show the legend
)

# Show the plot
fig.show()