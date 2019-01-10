def check_in_array(x, ref_value, fragment_threshold=5):
    try:
        return np.min(abs(x[:,0] - ref_value))/ref_value < fragment_threshold
    except:
        return np.nan