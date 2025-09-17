import math
import numpy as np
import sys
import pandas as pd
from scipy.spatial import ConvexHull
import random

np.set_printoptions(threshold=sys.maxsize)

'''
A collection of helper functions to extract and manipulate pixel information from the image
'''

#gets all the connected black pixels
def get_connected_black_pixels(image_np, start_x_np, start_y_np):

    '''
    The x and y values of the image must be switched because numpy calls its rows and columns differently,
    and they must be switched back in the end to match the CellProfiler coordinates.
    '''

    if image_np[start_x_np,start_y_np] > 0:
        return None

    else:
        x_img = len(image_np)
        y_img = len(image_np[0])

        directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

        stack = [(start_x_np, start_y_np)]

        visited = set()
        visited.add((start_x_np, start_y_np))

        connected_pixels = [(start_x_np, start_y_np)]

        while stack:
            x, y = stack.pop()

            for dx, dy in directions:
                nx, ny = x + dx, y + dy

                if 0 <= nx < (x_img - 1) and 0 <= ny < (y_img - 1) and image_np[nx,ny] == 0:

                    if (nx, ny) not in visited:
                        visited.add((nx, ny))
                        stack.append((nx, ny))
                        connected_pixels.append((nx, ny))

        connected_pixels_img = [(pixel_position[1], pixel_position[0]) for pixel_position in connected_pixels]

        return(connected_pixels_img)

#finds the eges of a set of connected pixels
def find_edge_pixels(connected_pixels, method = 'edge'):

    if method == 'edge':

        edge_pixels = []
        edge_neighbouring_pixels = set()
        for pixel in connected_pixels:
            neighbours = [(-1, 0), (1, 0), (0, -1), (0, 1), (1, 1), (-1, 1), (1, -1), (-1, -1)]
            neighbouring_pixels = []

            for directions in neighbours:
                neighbouring_pixels.append((pixel[0] + directions[0], pixel[1] + directions[1]))

            if all(element in connected_pixels for element in neighbouring_pixels) == False:
                edge_pixels.append(pixel)

                for element in list(set(neighbouring_pixels).difference(connected_pixels)):
                    edge_neighbouring_pixels.add(element)

    elif method == 'ConvexHull':
        connected_pixels_np = np.asarray(connected_pixels)
        hull = ConvexHull(connected_pixels_np)

        edge_pixels = connected_pixels_np[hull.vertices]

    elif method == 'ConcaveHull':
        #https://gist.github.com/jclosure/d93f39a6c7b1f24f8b92252800182889
        print('lets see')

    return(list(edge_neighbouring_pixels))

#organizes the pixels so that they can be made into a polygon
def pixel_sort(pixels_to_sort, method = 'frog_jump'):

    try:

        if method == 'frog_jump':

            directions_pixel_one = [(-1, 0), (1, 0), (0, -1), (0, 1)]
            coors_sorted = []
            next_pixel = []

            while len(coors_sorted) == 0:

                coor = random.choice(pixels_to_sort)
                neighbouring_pixels = []

                for direction in directions_pixel_one:
                    neighbouring_pixels.append((coor[0] + direction[0], coor[1] + direction[1]))

                coor_counts = 0
                for neighbouring_pixel in neighbouring_pixels:

                    if neighbouring_pixel in pixels_to_sort:
                        coor_counts = coor_counts + 1

                if coor_counts == 2:
                    for neighbouring_pixel in neighbouring_pixels:

                        if neighbouring_pixel in pixels_to_sort:
                            next_pixel.append(neighbouring_pixel)
                            coors_sorted.append(coor)
                            break

            while next_pixel:
                start_pixel = next_pixel.pop()

                if start_pixel[0] - coors_sorted[-1][0] == 1:
                    potential_next_neighbouring_pixel = (start_pixel[0] + 1, start_pixel[1])

                    if potential_next_neighbouring_pixel in pixels_to_sort and potential_next_neighbouring_pixel not in coors_sorted:
                        coors_sorted.append(start_pixel)
                        next_pixel.append(potential_next_neighbouring_pixel)

                    elif potential_next_neighbouring_pixel in pixels_to_sort and potential_next_neighbouring_pixel in coors_sorted:
                        coors_sorted.append(start_pixel)

                    else:
                        side_directions = [(0, -1), (0, 1)]
                        new_directions = []

                        for direction in side_directions:
                            new_directions.append((start_pixel[0] + direction[0], start_pixel[1] + direction[1]))

                        for new_direction in new_directions:

                            if new_direction in pixels_to_sort and new_direction not in coors_sorted:
                                coors_sorted.append(start_pixel)
                                next_pixel.append(new_direction)

                            elif new_direction in pixels_to_sort and new_direction in coors_sorted:
                                coors_sorted.append(start_pixel)

                if start_pixel[0] - coors_sorted[-1][0] == -1:
                    potential_next_neighbouring_pixel = (start_pixel[0] - 1, start_pixel[1])

                    if potential_next_neighbouring_pixel in pixels_to_sort and potential_next_neighbouring_pixel not in coors_sorted:
                        coors_sorted.append(start_pixel)
                        next_pixel.append(potential_next_neighbouring_pixel)

                    elif potential_next_neighbouring_pixel in pixels_to_sort and potential_next_neighbouring_pixel in coors_sorted:
                        coors_sorted.append(start_pixel)

                    else:
                        side_directions = [(0, -1), (0, 1)]
                        new_directions = []

                        for direction in side_directions:
                            new_directions.append((start_pixel[0] + direction[0], start_pixel[1] + direction[1]))

                        for new_direction in new_directions:

                            if new_direction in pixels_to_sort and new_direction not in coors_sorted:
                                coors_sorted.append(start_pixel)
                                next_pixel.append(new_direction)

                            elif new_direction in pixels_to_sort and new_direction in coors_sorted:
                                coors_sorted.append(start_pixel)

                if start_pixel[1] - coors_sorted[-1][1] == 1:
                    potential_next_neighbouring_pixel = (start_pixel[0], start_pixel[1] + 1)

                    if potential_next_neighbouring_pixel in pixels_to_sort and potential_next_neighbouring_pixel not in coors_sorted:
                        coors_sorted.append(start_pixel)
                        next_pixel.append(potential_next_neighbouring_pixel)

                    elif potential_next_neighbouring_pixel in pixels_to_sort and potential_next_neighbouring_pixel in coors_sorted:
                        coors_sorted.append(start_pixel)

                    else:
                        side_directions = [(-1, 0), (1, 0)]
                        new_directions = []

                        for direction in side_directions:
                            new_directions.append((start_pixel[0] + direction[0], start_pixel[1] + direction[1]))

                        for new_direction in new_directions:

                            if new_direction in pixels_to_sort and new_direction not in coors_sorted:
                                coors_sorted.append(start_pixel)
                                next_pixel.append(new_direction)

                            elif new_direction in pixels_to_sort and new_direction in coors_sorted:
                                coors_sorted.append(start_pixel)

                if start_pixel[1] - coors_sorted[-1][1] == -1:
                    potential_next_neighbouring_pixel = (start_pixel[0], start_pixel[1] - 1)

                    if potential_next_neighbouring_pixel in pixels_to_sort and potential_next_neighbouring_pixel not in coors_sorted:
                        coors_sorted.append(start_pixel)
                        next_pixel.append(potential_next_neighbouring_pixel)

                    elif potential_next_neighbouring_pixel in pixels_to_sort and potential_next_neighbouring_pixel in coors_sorted:
                        coors_sorted.append(start_pixel)

                    else:
                        side_directions = [(-1, 0), (1, 0)]
                        new_directions = []

                        for direction in side_directions:
                            new_directions.append((start_pixel[0] + direction[0], start_pixel[1] + direction[1]))

                        for new_direction in new_directions:

                            if new_direction in pixels_to_sort and new_direction not in coors_sorted:
                                coors_sorted.append(start_pixel)
                                next_pixel.append(new_direction)

                            elif new_direction in pixels_to_sort and new_direction in coors_sorted:
                                coors_sorted.append(start_pixel)

        elif method == 'clockwise':

            centre_x, centre_y = sum([x for x, _ in pixels_to_sort]) / len(pixels_to_sort), sum([y for _, y in pixels_to_sort]) / len(pixels_to_sort)

            angles = [math.atan2(y - centre_y, x - centre_x) for x, y in pixels_to_sort]

            clockwise_indices = sorted(range(len(pixels_to_sort)), key=lambda i: angles[i])
            coors_sorted = [pixels_to_sort[i] for i in clockwise_indices]

        if len(coors_sorted) >= len(pixels_to_sort) * 0.75:
            return (coors_sorted)

        else:
            return pixel_sort(pixels_to_sort)

    except RecursionError:
        
        return(None)