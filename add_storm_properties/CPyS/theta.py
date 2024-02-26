import numpy as np

def theta(x0=120, x1=130, y0=12, y1=10):
    """
    Computes the angular direction between two points.
    0° corresponds to eastward, 90° northward, 180° westward and 270° southward.

    Parameters
    ----------
    x0: longitude coordinate of the current point
    x1: longitude coordinate of the next point
    y0: latitude coordinate of the current point
    y1: longitude coordinate of the next point

    Returns
    -------
    th (float): The directional angle between the current and the next point, in degrees.
    """
    v = [x1 - x0, y1 - y0] # Vector corresponding to the track movement
    if v != [0,0]: # If not stagnating
        cos = (x1 - x0) / ( np.linalg.norm(v) )# Cosinus with eastward vector [1,0]
        if cos == -1:
            th = 180
        else:
            th = np.sign(y1 - y0) * np.arccos(cos) * 180 / np.pi
    else: # If stagnating: Set to NaN
        th = np.nan
    if th < 0: th += 360; # Make all angles between 0 and 360
    return th


def theta_track(lon, lat):
    """
    Computes the angular direction for each points along a track.
    Handling the track altogether allows for treating stationnary cases as well as the end of the track.
    In stationnary cases and in the end of the track, the direction from the previous point is taken for the current point.

    Parameters
    ----------
    lon: The list of longitudes along the track
    lat: The list of latitude along the track

    Returns
    -------
    th (list): values of th along the track.
    """
    th = []
    assert len(lon) == len(lat), "The two vector do not have the same length"
    n = len(lon)

    # Compute the angle for all the points
    for i in range(
            n - 1
    ):  # Computing the direction between each point and the following
        th.append(theta(lon[i], lon[i + 1], lat[i], lat[i + 1]))
        if np.isnan(th[-1]) & (
                i != 0
        ):  # If two successive points are superimposed, we take the previous direction
            th[-1] = th[-2]
    if n > 1:
        th.append(th[-1])
        # The direction for the last point is considered the same as the point before
    else:
        th = [np.nan]
    return th


def theta_multitrack(tracks):
    """
    Compute the angular direction for every tracks in a dataset.
    All tracks must have at least two points.

    Parameters
    ----------
    tracks (pd.DataFrame): The set of TC points including columns:
        * time
        * lon
        * lat

    Returns
    -------
    thetas (list): The list of angle for each point in the dataset
    """

    assert tracks.groupby("track_id").time.count().min() > 1, "The dataset contains tracks with only one point."

    # Identify position of track point from the end of the track (tpos=0 for the last point)
    tracks["tpos"] = tracks.sort_values("time", ascending=False).groupby("track_id").transform("cumcount")

    lon, lat = tracks.lon.values, tracks.lat.values

    # Run theta_track as if it was all one track
    th = np.array(theta_track(lon, lat))

    # Manage last point of each track : Set same angle as the point before
    th[list((tracks.tpos == 0).values)] = th[list((tracks.tpos == 1).values)]

    return th

if __name__ == "__main__":
    import numpy as np
    import pandas as pd

    # Test theta
    assert theta(0,1,0,0) == 0.0 # Eastward
    assert theta(0,0,0,1) == 90.0 # Northward
    assert theta(0,-1,0,0) == 180.0 # Westward
    assert theta(0,0,0,-1) == 270.0 # Southward

    # Test theta_track
    assert theta_track([0,1,1,0,0], [0,0,1,1,0]) == [0.0, 90.0, 180, 270.0, 270.0]

    # Test theta_multitrack
    tracks = pd.read_csv("tests/1996.csv", index_col = False)
    assert type(theta_multitrack(tracks)) == np.ndarray
    assert theta_multitrack(tracks)[1] == 90.

    print("All good")
