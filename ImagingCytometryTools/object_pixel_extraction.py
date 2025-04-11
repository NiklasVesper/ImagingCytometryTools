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
def get_connected_black_pixels(image_np, start_x_np, start_y_np): #here the x and y values of the image must be switched because numpy calls its rows and columns differently

    #checks if the starting pixel is white, so one does not jump outside the cell boundaries
    if image_np[start_x_np,start_y_np] > 0:
        return None

    else:
        x_img = len(image_np) #x dimensions for the image
        y_img = len(image_np[0]) #y dimensions for the image

        directions = [(-1, 0), (1, 0), (0, -1), (0, 1)] #serch directions ((1,1) etc. are not an option because it would jump over the cell bondaries)

        stack = [(start_x_np, start_y_np)] #stack list for the seach

        visited = set() #set for seched pixels
        visited.add((start_x_np, start_y_np))

        connected_pixels = [(start_x_np, start_y_np)] # list for all the connectec pixels

        # while there is a stack run the following
        while stack:
            x, y = stack.pop() #removes the last element of the stack

            #goes throw all the elements in the directions
            for dx, dy in directions:
                nx, ny = x + dx, y + dy #creates a pair of new pixel positions

                if 0 <= nx < (x_img - 1) and 0 <= ny < (y_img - 1) and image_np[nx,ny] == 0: #checks if the new position is a black pixel and in the image (just to be shure)

                    if (nx, ny) not in visited: #checks if the new position has be visited
                        visited.add((nx, ny)) #appends to visited
                        stack.append((nx, ny)) #appends to the stack
                        connected_pixels.append((nx, ny)) #appends to the connected pixels

        connected_pixels_img = [(pixel_position[1], pixel_position[0]) for pixel_position in connected_pixels] #switches the x and y position to match the CellProfiler coordinates

        return(connected_pixels_img)

#finds the eges of a set of connected pixels
def find_edge_pixels(connected_pixels, method = 'edge'):

    if method == 'edge': #does not ignore wholes in the connected area

        edge_pixels = [] #empty list for the edge pixels
        edge_neighbouring_pixels = set()
        for pixel in connected_pixels:
            neighbours = [(-1, 0), (1, 0), (0, -1), (0, 1), (1, 1), (-1, 1), (1, -1), (-1, -1)] #movement list for the neighboring pixels
            neighbouring_pixels = [] #empty list for the pixels

            #creates all the potential neighboring pixels to be checked
            for directions in neighbours:
                neighbouring_pixels.append((pixel[0] + directions[0], pixel[1] + directions[1]))

            #initializes the search condition so not all have to be checked
            if all(element in connected_pixels for element in neighbouring_pixels) == False:
                edge_pixels.append(pixel)

                #checks if the potential neighboring pixels are actual neighboring pixels
                for element in list(set(neighbouring_pixels).difference(connected_pixels)):
                    edge_neighbouring_pixels.add(element)

    elif method == 'ConvexHull': #does ignore wholes in the connected area
        connected_pixels_np = np.asarray(connected_pixels) #turns the connected pixels into a numpy array
        hull = ConvexHull(connected_pixels_np) #calculates convex hull

        edge_pixels = connected_pixels_np[hull.vertices]

    elif method == 'ConcaveHull':
        #https://gist.github.com/jclosure/d93f39a6c7b1f24f8b92252800182889
        print('')

    return(list(edge_neighbouring_pixels))

#organizes the pixels so that they can be made into a polygon
def pixel_sort(pixels_to_sort, method = 'frog_jump'):

    #just jumps from one direction to another memorizing the last direction it came form
    if method == 'frog_jump':

        directions_pixel_one = [(-1, 0), (1, 0), (0, -1), (0, 1)] #allowed jump directions
        coors_sorted = [] #empty list for the sorted coordinates
        next_pixel = [] #empty list for the next pixel

        #initializes the search if there are no coordinates in the coors_sorted list
        while len(coors_sorted) == 0:

            coor = random.choice(pixels_to_sort) #picks a random start coordinate
            neighbouring_pixels = [] #empty list for the neighbouring pixels

            #creates the neighboring pixels
            for direction in directions_pixel_one:
                neighbouring_pixels.append((coor[0] + direction[0], coor[1] + direction[1]))

            #checks whether only two neighboring coordinates are found viua a counter
            coor_counts = 0
            for neighbouring_pixel in neighbouring_pixels:

                if neighbouring_pixel in pixels_to_sort:
                    coor_counts = coor_counts + 1

            #if there are two found the while loop is broken
            if coor_counts == 2:
                for neighbouring_pixel in neighbouring_pixels:

                    #the found coordinates are appended into the list to continue the search
                    if neighbouring_pixel in pixels_to_sort:
                        next_pixel.append(neighbouring_pixel)
                        coors_sorted.append(coor)
                        break

        #runs while there is a next pixel
        while next_pixel:
            start_pixel = next_pixel.pop()

            if start_pixel[0] - coors_sorted[-1][0] == 1: #gets the direction via the pixel difference
                potential_next_neighbouring_pixel = (start_pixel[0] + 1, start_pixel[1]) #generates the potential next pixel

                #checks if the potential next pixel exists within the list
                if potential_next_neighbouring_pixel in pixels_to_sort and potential_next_neighbouring_pixel not in coors_sorted:
                    coors_sorted.append(start_pixel)
                    next_pixel.append(potential_next_neighbouring_pixel)

                #checks whether the pixel is already in the sorted list
                elif potential_next_neighbouring_pixel in pixels_to_sort and potential_next_neighbouring_pixel in coors_sorted:
                    coors_sorted.append(start_pixel)

                #generates the side coordinates to check on other directions if no pixel infrot is found
                else:
                    side_directions = [(0, -1), (0, 1)] #directions to go to
                    new_directions = [] #empty list for the new directions

                    #generates the directions
                    for direction in side_directions:
                        new_directions.append((start_pixel[0] + direction[0], start_pixel[1] + direction[1]))

                    #goes throw the directions
                    for new_direction in new_directions:

                        # checks if the potential next pixel exists within the list
                        if new_direction in pixels_to_sort and new_direction not in coors_sorted:
                            coors_sorted.append(start_pixel)
                            next_pixel.append(new_direction)

                        # checks whether the pixel is already in the sorted list
                        elif new_direction in pixels_to_sort and new_direction in coors_sorted:
                            coors_sorted.append(start_pixel)

            if start_pixel[0] - coors_sorted[-1][0] == -1: #gets the direction via the pixel difference
                potential_next_neighbouring_pixel = (start_pixel[0] - 1, start_pixel[1]) #generates the potential next pixel

                # checks if the potential next pixel exists within the list
                if potential_next_neighbouring_pixel in pixels_to_sort and potential_next_neighbouring_pixel not in coors_sorted:
                    coors_sorted.append(start_pixel)
                    next_pixel.append(potential_next_neighbouring_pixel)

                # checks whether the pixel is already in the sorted list
                elif potential_next_neighbouring_pixel in pixels_to_sort and potential_next_neighbouring_pixel in coors_sorted:
                    coors_sorted.append(start_pixel)

                # generates the side coordinates to check on other directions if no pixel infrot is found
                else:
                    side_directions = [(0, -1), (0, 1)] #directions to go to
                    new_directions = [] #empty list for the new directions

                    # generates the directions
                    for direction in side_directions:
                        new_directions.append((start_pixel[0] + direction[0], start_pixel[1] + direction[1]))

                    # goes throw the directions
                    for new_direction in new_directions:

                        # checks if the potential next pixel exists within the list
                        if new_direction in pixels_to_sort and new_direction not in coors_sorted:
                            coors_sorted.append(start_pixel)
                            next_pixel.append(new_direction)

                        # checks whether the pixel is already in the sorted list
                        elif new_direction in pixels_to_sort and new_direction in coors_sorted:
                            coors_sorted.append(start_pixel)

            if start_pixel[1] - coors_sorted[-1][1] == 1: #gets the direction via the pixel difference
                potential_next_neighbouring_pixel = (start_pixel[0], start_pixel[1] + 1) #generates the potential next pixel

                # checks if the potential next pixel exists within the list
                if potential_next_neighbouring_pixel in pixels_to_sort and potential_next_neighbouring_pixel not in coors_sorted:
                    coors_sorted.append(start_pixel)
                    next_pixel.append(potential_next_neighbouring_pixel)

                elif potential_next_neighbouring_pixel in pixels_to_sort and potential_next_neighbouring_pixel in coors_sorted:
                    coors_sorted.append(start_pixel)

                # generates the side coordinates to check on other directions if no pixel infrot is found
                else:
                    side_directions = [(-1, 0), (1, 0)] #directions to go to
                    new_directions = [] #empty list for the new directions

                    # generates the directions
                    for direction in side_directions:
                        new_directions.append((start_pixel[0] + direction[0], start_pixel[1] + direction[1]))

                    # goes throw the directions
                    for new_direction in new_directions:

                        # checks if the potential next pixel exists within the list
                        if new_direction in pixels_to_sort and new_direction not in coors_sorted:
                            coors_sorted.append(start_pixel)
                            next_pixel.append(new_direction)

                        # checks whether the pixel is already in the sorted list
                        elif new_direction in pixels_to_sort and new_direction in coors_sorted:
                            coors_sorted.append(start_pixel)

            if start_pixel[1] - coors_sorted[-1][1] == -1: #gets the direction via the pixel difference
                potential_next_neighbouring_pixel = (start_pixel[0], start_pixel[1] - 1) #generates the potential next pixel

                # checks if the potential next pixel exists within the list
                if potential_next_neighbouring_pixel in pixels_to_sort and potential_next_neighbouring_pixel not in coors_sorted:
                    coors_sorted.append(start_pixel)
                    next_pixel.append(potential_next_neighbouring_pixel)

                elif potential_next_neighbouring_pixel in pixels_to_sort and potential_next_neighbouring_pixel in coors_sorted:
                    coors_sorted.append(start_pixel)

                # generates the side coordinates to check on other directions if no pixel infrot is found
                else:
                    side_directions = [(-1, 0), (1, 0)] #directions to go to
                    new_directions = [] #empty list for the new directions

                    # generates the directions
                    for direction in side_directions:
                        new_directions.append((start_pixel[0] + direction[0], start_pixel[1] + direction[1]))

                    # goes throw the directions
                    for new_direction in new_directions:

                        # checks if the potential next pixel exists within the list
                        if new_direction in pixels_to_sort and new_direction not in coors_sorted:
                            coors_sorted.append(start_pixel)
                            next_pixel.append(new_direction)

                        # checks whether the pixel is already in the sorted list
                        elif new_direction in pixels_to_sort and new_direction in coors_sorted:
                            coors_sorted.append(start_pixel)

    #sorts the pixels clockwise around a centroid
    elif method == 'clockwise':

        centre_x, centre_y = sum([x for x, _ in pixels_to_sort]) / len(pixels_to_sort), sum([y for _, y in pixels_to_sort]) / len(pixels_to_sort)  # calculates center of the area
    
        angles = [math.atan2(y - centre_y, x - centre_x) for x, y in pixels_to_sort]  # gives the angels between the center and the pixel positions
    
        #gives the indices of the positions and sorts the pixels
        clockwise_indices = sorted(range(len(pixels_to_sort)), key=lambda i: angles[i])
        coors_sorted = [pixels_to_sort[i] for i in clockwise_indices]

    #returns the sorted pixels if the correct amount of pixels is found
    if len(coors_sorted) >= len(pixels_to_sort) * 0.75: # maybe bug sorce
        return (coors_sorted)

    #reruns the function if not the correct amount of pixels is found
    else:
        return pixel_sort(pixels_to_sort)