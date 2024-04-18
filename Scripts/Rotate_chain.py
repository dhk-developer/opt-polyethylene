import numpy as np

def Rotator(AP_chain):
    """
    Randomly rotates a Polyethylene chain to increase variation and non-uniformity of initial structure.

    :param AP_chain: Single PE chain within box.

    :return: Randomly rotated Polyethylene chain - 0 to 360 degrees.
    """

    n_AP_chain = []

    # random angle
    rand_angle = np.random.uniform(0, 2*np.pi)

    print ("works so far")

    # Compute cosine and sine of the random angle
    cos_angle = np.cos(rand_angle)
    sin_angle = np.sin(rand_angle)

    # Define the rotation matrix
    R_yz = np.array([[1, 0, 0],
                     [0, cos_angle, sin_angle],
                     [0, -sin_angle, cos_angle]])


    # For each atom's atomic position in the chain, apply the rotation matrix to the coordinates
    for i in AP_chain:
        print(i)
        n_AP_chain.append(np.dot(i, R_yz))

    return n_AP_chain




