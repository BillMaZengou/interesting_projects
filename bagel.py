import numpy as np
from numpy import cos as C, sin as S
import os

class Rotation(object):
    """
    Give the axis 'x', 'y', or 'z', return a rotation matrix
    """
    def __init__(self, axis):
        super(Rotation, self).__init__()
        self.axis = axis  # rotation around this axis

    def rotation_matrix(self, angle):
        if self.axis == 'x':
            self.r_matrix = np.array(
                [
                    [1,         0,         0],
                    [0,  C(angle), -S(angle)],
                    [0,  S(angle),  C(angle)]
                                                ]
            )
        elif self.axis == 'y':
            self.r_matrix = np.array(
                [
                    [ C(angle), 0, S(angle)],
                    [        0, 1,        0],
                    [-S(angle), 0, C(angle)]
                                                ]
            )
        elif self.axis == 'z':
            self.r_matrix = np.array(
                [
                    [ C(angle), -S(angle), 0],
                    [ S(angle),  C(angle), 0],
                    [        0,         0, 1]
                                                ]
            )
        else:
            self.r_matrix = np.zeros((3, 3))
            print("ERROR: only x, y, z allowed")
        return self.r_matrix

class Point_Light(object):
    def __init__(self, x, y, z):
        super(Point_Light, self).__init__()
        self.position = np.array([x, y, z])
        self.position.shape = (3, 1)

    def origin_to_light(self):
        return self.position

class Donut_Point(object):
    def __init__(self, R1, R2, theta, phi, distance):
        super(Donut_Point, self).__init__()
        self.theta = theta
        self.phi = phi

        # Construct the donut
        self.centre = np.array([R2, 0, 0])
        self.centre.shape = (3, 1)
        self.circle = np.array([R1*C(self.theta), R1*S(self.theta), 0])
        self.circle.shape = (3, 1)
        self.revolution = Rotation('y')
        self.translation = np.array([0, 0, distance])
        self.translation.shape = (3, 1)
        self.donut = self.revolution.rotation_matrix(self.phi) @ (self.centre + self.circle) + self.translation

        # Transformation for visual apperence
        self.x_rotation = Rotation('x')
        self.z_rotation = Rotation('z')

        # Normal Line
        self.normal_direction = np.array([C(self.theta), S(self.theta), 0])
        self.normal_direction.shape = (3, 1)
        self.normal_direction = self.revolution.rotation_matrix(self.phi) @ self.normal_direction + self.translation

    def rotate_donut(self, z_angle, x_angle):
        return self.x_rotation.rotation_matrix(x_angle) @ self.z_rotation.rotation_matrix(z_angle) @ (self.donut - self.translation) + self.translation

    def rotate_normal(self, z_angle, x_angle):
        return self.x_rotation.rotation_matrix(x_angle) @ self.z_rotation.rotation_matrix(z_angle) @ (self.normal_direction - self.translation) + self.translation

def main():
    """----- Camera -------"""
    aspect_ratio = 16 / 9
    screen_height = 40
    screen_width = int(aspect_ratio * screen_height)

    output = []
    z_buffer = []
    for j in range(screen_height):
        output.append([])
        z_buffer.append([])
        for i in range(screen_width):
            output[j].append(' ')
            z_buffer[j].append(0)

    """----- Light -------"""
    point_light = Point_Light(0, 1, -1)

    """----- Model -------"""
    R1 = 2                  # cross-sectional radius
    R2 = 6                  # torus radius
    theta_spacing = 0.3     # cross-sectional angle
    phi_spacing = 0.06      # revolution angle
    donut_pos = []
    intensities = []

    # origin to model
    obj_distance = 20
    # origin to image plane
    z_prime = screen_height*obj_distance*6/(8*2*(R1+R2)) if screen_height < screen_width else screen_width*obj_distance*6/(8*2*(R1+R2))

    z_angle = 0
    x_angle = np.pi/4
    for theta in np.arange(0, 2*np.pi, theta_spacing):
        for phi in np.arange(0, 2*np.pi, phi_spacing):
            donut = Donut_Point(R1, R2, theta, phi, obj_distance)
            donut_point = donut.rotate_donut(z_angle, x_angle)
            donut_normal = donut.rotate_normal(z_angle, x_angle)
            intensities.append(donut_normal.T @ point_light.origin_to_light())
            donut_pos.append(donut_point)
    num_points = len(intensities)

    """----- Update -------"""
    while True:
        # Rescale to 0 - 11
        intensities = (intensities - min(intensities))/(max(intensities) - min(intensities)) * 11

        for i in range(num_points):
            donut_point = donut_pos[i]
            x = donut_point[0][0]
            y = donut_point[1][0]
            reciprocal_z = 1 / donut_point[2][0]
            x_pixel = int(z_prime*x*reciprocal_z*aspect_ratio + screen_width/2)
            y_pixel = int(-z_prime*y*reciprocal_z + screen_height/2)
            l = int(intensities[i])
            if (x_pixel > 0 and x_pixel < screen_width and y_pixel > 0 and y_pixel < screen_height):
                if reciprocal_z > z_buffer[y_pixel][x_pixel]:
                    output[y_pixel][x_pixel] = ".,-~:;=!*#$@"[l]
                    z_buffer[y_pixel][x_pixel] = reciprocal_z

        """----- Display -------"""
        for j in range(len(output)):
            for i in range(len(output[0])):
                if i < (screen_width-1):
                    print(output[j][i], end='')
                else:
                    print(output[j][i])

        """----- Next Frame -------"""
        z_angle += 0.06
        x_angle += 0.02

        temp = 0
        for theta in np.arange(0, 2*np.pi, theta_spacing):
            for phi in np.arange(0, 2*np.pi, phi_spacing):
                donut = Donut_Point(R1, R2, theta, phi, obj_distance)
                donut_point = donut.rotate_donut(z_angle, x_angle)
                donut_normal = donut.rotate_normal(z_angle, x_angle)
                intensities[temp] = donut_normal.T @ point_light.origin_to_light()
                donut_pos[temp] = donut_point
                temp += 1

        """----- Clear -------"""
        for j in range(screen_height):
            for i in range(screen_width):
                output[j][i] = ' '
                z_buffer[j][i] = 0
        _ = os.system("clear")


if __name__ == '__main__':
    main()
