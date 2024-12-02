import dash
from dash import dcc, html
from dash.dependencies import Input, Output
import plotly.graph_objs as go
import serial
import datetime

# Initialize the serial port
# Replace 'COM3' with your actual port (e.g., '/dev/ttyUSB0' on Linux, 'COM3' on Windows)
SERIAL_PORT = 'COM8'
BAUD_RATE = 9600

try:
    ser = serial.Serial(SERIAL_PORT, BAUD_RATE, timeout=1)
except Exception as e:
    print(f"Error opening serial port: {e}")
    exit()

# Initialize the Dash app
app = dash.Dash(__name__, external_stylesheets=['https://codepen.io/chriddyp/pen/bWLwgP.css'])
app.title = "Real-Time Data Visualization from Serial Port"

# Initialize a list to store the data
times = []
values = []

# Layout of the dashboard
app.layout = html.Div([
    html.H1("Real-Time Data Visualization from Serial Port", style={'text-align': 'center'}),
    dcc.Graph(id='live-graph', style={'height': '65vh'}),
    dcc.Interval(
        id='interval-component',
        interval=1000,  # Update every 1000 ms (1 second)
        n_intervals=0  # Counter for intervals
    )
])

# Update the plot every interval
@app.callback(
    Output('live-graph', 'figure'),
    Input('interval-component', 'n_intervals')
)
def update_graph(n):
    global ser

    # Read data from the serial port
    try:
        if ser.in_waiting > 0:
            line = ser.readline().decode('utf-8').strip()  # Read a line of data
            new_value = float(line)  # Convert it to a float (adjust according to your data format)
        else:
            new_value = None
    except Exception as e:
        print(f"Error reading serial data: {e}")
        new_value = None

    # Update the graph only if new data is received
    if new_value is not None:
        current_time = datetime.datetime.now().strftime('%H:%M:%S')
        times.append(current_time)
        values.append(new_value)

        # Keep the last 100 data points for display
        if len(times) > 100:
            times.pop(0)
            values.pop(0)

    # Create a Plotly graph
    fig = go.Figure(
        data=[
            go.Scatter(
                x=times,
                y=values,
                mode='lines+markers',
                line=dict(color='royalblue', width=2),
                marker=dict(size=5),
                name='Serial Data'
            )
        ]
    )

    # Customize the layout
    fig.update_layout(
        title="Real-Time Data from Serial Port",
        xaxis_title="Time",
        yaxis_title="Value",
        xaxis=dict(showline=True, showgrid=True),
        yaxis=dict(showline=True, showgrid=True),
        template='plotly_dark'
    )

    return fig

# Run the Dash app
if __name__ == '__main__':
    app.run_server(debug=True)
