#-- my_code_hw01.py
#-- hw01 GEO1015.2022
#-- [YOUR NAME]
#-- [YOUR STUDENT NUMBER]
import random
import numpy as np
import math
from scipy.spatial import *
import matplotlib.pyplot as plt
import startinpy
import time


def nn_xy(dt, kd, all_z, x, y):
    """
    !!! TO BE COMPLETED !!!

    Function that interpolates with nearest neighbour method.

    Input:
        dt:     the DT of the input points (a startinpy object)
        kd:     the kd-tree of the input points
        all_z:  an array with all the z values, same order as kd.data
        x:      x-coordinate of the interpolation location
        y:      y-coordinate of the interpolation location
    Output:
        z: the estimation of the height value,
           (raise Exception if outside convex hull)
    """
    # the task is just to find the nearest point
    # -- kd-tree docs: https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.html#scipy.spatial.KDTree
    # -- you are *not* allowed to use the function for the nn interpolation that I wrote for startinpy
    # -- you need to write your own code for this step
    # z = random.uniform(0, 100)
    # raise Exception("Outside convex hull")
    try:
        vertex = dt.closest_point(x, y)
        getted_pt = dt.get_point(vertex)
        z = getted_pt[-1]
    except Exception as e:
        print(e)
        raise Exception("Outside convex hull")
    z_hugo = dt.interpolate_nn(x,y)
    #print("z_hugo and z_max is",Z,z_hugo,x,y)
    if abs(z-z_hugo) > 0.001:
        print("an erro is ",abs(z-z_hugo))

    return z


def idw_xy(dt, kd, all_z, x, y, power, radius):
    """
    !!! TO BE COMPLETED !!!

    Function that interpolates with IDW

    Input:
        dt:     the DT of the input points (a startinpy object)
        kd:     the kd-tree of the input points
        all_z:  an array with all the z values, same order as kd.data
        x:      x-coordinate of the interpolation location
        y:      y-coordinate of the interpolation location
        power:  power to use for IDW
        radius: search radius
¨    Output:
        z: the estimation of the height value,
           (raise Exception if (1) outside convex hull or (2) no point in search radius
    """
    # -- kd-tree docs: https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.html#scipy.spatial.KDTree
    if dt.is_inside_convex_hull(x, y) != True:
        raise Exception("Outside convex boundary")
    pts = kd.query_ball_point([x, y], r=radius)
    #print(pts)
    nearest = []
    for i in pts:
        pt = kd.data[i]
        distance = math.sqrt((x - pt[0]) ** 2 + (y - pt[1]) ** 2)
        #print(dt.closest_point(*pt))
        z = all_z[dt.closest_point(*pt)-1]
        nearest.append((z, distance ** (-power)))
    if nearest == []:
        # print(x,y)
        raise Exception("No points in search radius")
    # print(nearest)
    values, weights = zip(*nearest)
    #print(values,weights)
    sums = sum(weights)
    # print(sums)
    z = 0

    for i in range(len(values)):
        z += (weights[i] / sums) * values[i]


    return z


def tin_xy(dt, kd, all_z, x, y):
    """
    !!! TO BE COMPLETED !!!

    Function that interpolates linearly in a TIN.

    Input:
        dt:     the DT of the input points (a startinpy object)
        kd:     the kd-tree of the input points
        all_z:  an array with all the z values, same order as kd.data
        x:      x-coordinate of the interpolation location
        y:      y-coordinate of the interpolation location
    Output:
        z: the estimation of the height value,
           (raise Exception if outside convex hull)
    """
    # -- startinpy docs: https://startinpy.rtfd.io/
    # -- you are *not* allowed to use the function for the TIN interpolation that I wrote for startinpy
    # -- you need to write your own code for this step
    try:
        tin = dt.locate(x, y)

        #print(tin)
    except:
        raise Exception("Outside convex hull")
    x1, y1, z1 = dt.get_point(tin[0])
    x2, y2, z2 = dt.get_point(tin[1])
    x3, y3, z3 = dt.get_point(tin[2])
    # AX=B A:[[X1,X2,X3],[Y1,Y2,Y3],[1,1,1]] X:[[W0],[W1],[W2]] B[[x],[y],[1]]
    X_matrix = np.array([[x1, x2, x3], [y1, y2, y3], [1, 1, 1]])
    Y_matrix = np.array([[x], [y], [1]])
    X_inv = np.linalg.inv(X_matrix)
    solution_matrix = np.dot(X_inv, Y_matrix)
    #print(solution_matrix)
    z = z1 * solution_matrix[0][0] + z2 * solution_matrix[1][0] + z3 * solution_matrix[2][0]
    # print(z)
    z_hugo = dt.interpolate_tin_linear(x,y)
    #print("z_hugo and z_max is",Z,z_hugo,x,y)
    if abs(z-z_hugo) > 0.001:
        print("an erro is ",abs(z-z_hugo))

    return z

def find_neighbor_points(Vor,pt_index,dt): # iterate ridge and put find vertices, neighbor point
    neibor_pair = Vor.ridge_points
    ridge_vertices = Vor.ridge_vertices
    #print("ridge vertices is",ridge_vertices)
    selected_point = []
    selected_vertices = []
    for i in range(len(neibor_pair)):
        if pt_index in neibor_pair[i]:
            selected_point.append(neibor_pair[i][neibor_pair[i] != pt_index][0])
            two_end_ridge = ridge_vertices[i]
            ava_pt = []
            for indice in two_end_ridge:
                pt = Vor.vertices[indice]
                if indice!=-1 and dt.is_inside_convex_hull(pt[0],pt[1]):
                    selected_vertices.append(pt)
                    ava_pt.append(True)
                else:
                    ava_pt.append(False)
            if True in ava_pt and False in ava_pt:

                #print("mid vertices",Vor.points[neibor_pair[i]])
                #print("end point is", Vor.vertices[two_end_ridge][ava_pt])
                two_pt = Vor.points[neibor_pair[i]]
                x1,y1 = 0.5*(two_pt[0]+two_pt[1])

                x2,y2 = Vor.vertices[two_end_ridge][ava_pt][0]
                try:

                    #points = modified_vertices_on_convex(dt,x1,y1,x2,y2)

                    modified_vertex_x,modified_vertex_y = modified_vertices_on_convex(dt,x1,y1,x2,y2)

                except Exception as e:

                    #hull = ConvexHull(dt.points[1:,:2])
                    voronoi_plot_2d(Vor)
                    #convex_hull_plot_2d(hull)
                    plt.scatter(x1, y1,c='g',s=100)
                    plt.scatter(x2, y2,c='r',s=100)
                    #plt.scatter(two_pt[0][0],two_pt[0][1],s = 100)
                    #plt.scatter(two_pt[1][0],two_pt[1][1],s = 100)
                    plt.show()
                    raise e
                selected_vertices.append(np.array([modified_vertex_x,modified_vertex_y]))
    #print("neighbor is ",selected_point)
    #print("vertices is ",selected_vertices)
    return selected_point,selected_vertices
def modified_vertices_on_convex(dt,inner_x,inner_y,outer_x,outer_y):
    convex = dt.convex_hull()
    x1 = inner_x-outer_x
    y1 = inner_y-outer_y

    vector1 = np.array([[x1],[y1]])
    vector0 = np.array([[outer_x], [outer_y]])
    points=[]
    segment= []
    for i in range(len(convex)-1):
        pt1 = dt.points[convex[i]][:2]
        pt2 = dt.points[convex[i+1]][:2]
        vector2 = pt2-pt1
        vector2 = vector2.reshape(2,1)
        vector3 = pt2.reshape(2,1)
        #print(vector1,vector2)
        vector_X = np.hstack([vector1,vector2])
        #print("x is",vector_X)
        vector_Y = vector3-vector0
        #points = []
        try:
            solution = np.linalg.inv(vector_X).dot(vector_Y)
            if 0<=solution[1][0] <= 1:
                #points.append(solution)
                #print("has solution",i)
                #print("vector3,vector2,solution[1][0]",vector3,vector2,solution[1][0],vector3-vector2*solution[1][0])
                point_x,point_y = vector3-vector2*solution[1][0] #shape 2,1
                #return point_x[0],point_y[0]
                points.append([point_x[0],point_y[0]])
                #segment.append(([pt1[0],pt2[0]],[pt1[1], pt2[1]]))
                #print(solution)
                #print(i)
        except Exception as e:
            print(e)
        #print(vector_X)
    #w1(x1,y1)+(inner_x,inner_y) = w2(x2,y2)+(x,y)
    #
    #print(type(points[0]))
    if len(points)==2:
        dist1 = (points[0][0]-outer_x)**2+(points[0][1]-outer_y)**2
        dist2 = (points[1][0] - outer_x) ** 2 + (points[1][1] - outer_y) ** 2
        if dist1>dist2:
            return points[0][0],points[0][1]
        else:
            return points[1][0],points[1][1]
    elif len(points)==1:
        return points[0][0],points[0][1]

    else:
        #len(points) < 10000:
        print(inner_x, inner_y,outer_x,outer_y)
        hull = ConvexHull(dt.points[1:, :2])
        convex_hull_plot_2d(hull)
        for X, Y in points:
            plt.scatter(X, Y, s=60, c='r')
        plt.scatter(inner_x, inner_y, s=60, c='y')

        plt.scatter(outer_x, outer_y, s=200, c='g')
        for X,Y in segment:
            plt.plot(X,Y)
        plt.show()
        print("points is",points)
        return points

# one hundred correct
def polygon_area(points):
    #3print(list(points))
    dt = startinpy.DT()
    dt.insert(points)
    #print(points,dt.points)
    #print(points)
    area = 0
    for tri in dt.triangles:

        #print(dt.points[tri[0]])
        x1,y1,z = dt.points[tri[0]]
        x2,y2,z = dt.points[tri[1]]
        x3,y3,z = dt.points[tri[2]]
        area_martrix = np.array([[x1,x2,x3],[y1,y2,y3],[1,1,1]])
        area += abs(np.linalg.det(area_martrix))
    if area == 0:
        return 0
        pts = dt.points
        trs = dt.triangles
        print(pts,trs)
        # -- plot
        plt.triplot(pts[:, 0], pts[:, 1], trs)
        # -- the vertex "0" shouldn't be plotted, so start at 1
        plt.plot(pts[1:, 0], pts[1:, 1], 'o')
        plt.show()
    return area

def nni_xy(dt, kd, all_z, x, y):
    if dt.is_inside_convex_hull(x,y) == False:
        raise Exception("Outside Boundary")
    #start = time.time()

    old_vor = Voronoi(kd.data,incremental=True)
    #end = time.time()
    #print(end-start)
    old_vertices = old_vor.vertices
    old_regions = old_vor.regions
    old_point_regioin = old_vor.point_region
    #print("start to calculate")
    old_vor.add_points(np.array([[x,y]]))
    old_vor.close()
    new_vertices = old_vor.vertices
    ridge_relation = old_vor.ridge_points
    vertex_added = len(old_vor.points) - 1
    polygons = []
    z_values = []
    for i in range(len(ridge_relation)):
        if vertex_added in ridge_relation[i]:
            polygon = []
            inter_pts = ridge_relation[i]
            inter_p = inter_pts[inter_pts != vertex_added][0]
            inter_p_coord = old_vor.points[inter_p]
            inter_p_z = dt.points[dt.closest_point(inter_p_coord[0],inter_p_coord[1])][2]
            z_values.append(inter_p_z)
            inter_p_vertices = old_regions[old_point_regioin[inter_p]]
            inter_p_vertices_coord = [old_vertices[vert] for vert in inter_p_vertices \
                                      if vert != -1 and old_vertices[vert] not in new_vertices]
            #print("interpolation point",inter_p_coord,inter_p_z)
            if inter_p_vertices_coord==[]:
                #print("null interpoltation point",inter_p,x,y)
                continue
            polygon.extend(inter_p_vertices_coord)
            #print("inner vertices",polygon)
            #inter_p_coord = old_vor.points[inter_p]
            #print(inter_p_coord)

            vertices = old_vor.ridge_vertices[i]
            #print("vertices",vertices)


            if -1 not in vertices:
                vertices1, vertices2 = vertices
                vertices1_coord,vertices2_coord=old_vor.vertices[[vertices1,vertices2]]
                polygon.append(vertices1_coord)
                polygon.append(vertices2_coord)
                #print(vertices1_coord,vertices1_coord)
            else:
                vertices1 = vertices[vertices!=-1]
                vertices_1_coord = old_vor.vertices[vertices1]
                polygon.append(vertices_1_coord)
                mid_x,mid_y = 0.5*(x+inter_p_coord[0]),0.5*(y+inter_p_coord[1])

                vertices_2_coord=modified_vertices_on_convex(dt, vertices_1_coord[0],vertices_1_coord[1],mid_x,mid_y)
                polygon.append(vertices_2_coord)
                #print("vert",vertices1)
            polygons.append(polygon)
            #print("polygon",polygon)

            #print(new_vertices[vertices1],new_vertices[vertices2])
            #polygon.extend([new_vertices[vertices1],new_vertices[vertices2]])
    total_area = 0
    sum = 0
    #print(z_values,"asdfasd",polygons)
    for z, poly in zip(z_values,polygons):
        points = [[i[0],i[1],z] for i in poly]
        #print(points)
        area = polygon_area(points)
        total_area+=area
        sum+=z*area
        #print(area)
    Z = sum/total_area
    #print("Z value is ",Z)
    #start = time.time()
    #print("second roudn",start-end)
    z_hugo = dt.interpolate_nni(x,y)
    #print("z_hugo and z_max is",Z,z_hugo,x,y)
    if abs(Z-z_hugo) > 0.001:
        print("an erro is ",abs(Z-z_hugo))

    return Z
    #vertices = collect_vertices_vor_xy(dt,x,y)
def nni_xy_modified(dt, kd, all_z, x, y):
    if dt.is_inside_convex_hull(x,y) == False:
        #return -9999999
        raise Exception("Outside Boundary")
    pts = find_surroundingPoint(dt,x,y)[:,:2]
    old_vor = Voronoi(pts, incremental=1)
    old_vertices = old_vor.vertices
    old_regions = old_vor.regions
    old_point_regioin = old_vor.point_region
    # print("start to calculate")
    #print(pts)
    old_vor.add_points(np.array([[x, y]]))
    old_vor.close()
    new_vertices = old_vor.vertices
    ridge_relation = old_vor.ridge_points
    vertex_added = len(old_vor.points) - 1
    polygons = []
    z_values = []
    for i in range(len(ridge_relation)):
        if vertex_added in ridge_relation[i]:
            polygon = []
            inter_pts = ridge_relation[i]
            inter_p = inter_pts[inter_pts != vertex_added][0]
            inter_p_coord = old_vor.points[inter_p]
            inter_p_z = dt.points[dt.closest_point(inter_p_coord[0], inter_p_coord[1])][2]
            z_values.append(inter_p_z)
            inter_p_vertices = old_regions[old_point_regioin[inter_p]]
            inter_p_vertices_coord = [old_vertices[vert] for vert in inter_p_vertices \
                                      if vert != -1 and old_vertices[vert] not in new_vertices]
            # print("interpolation point",inter_p_coord,inter_p_z)
            if inter_p_vertices_coord == []:
                # print("null interpoltation point",inter_p,x,y)
                continue
            polygon.extend(inter_p_vertices_coord)
            # print("inner vertices",polygon)
            # inter_p_coord = old_vor.points[inter_p]
            # print(inter_p_coord)

            vertices = old_vor.ridge_vertices[i]
            # print("vertices",vertices)

            if -1 not in vertices:
                vertices1, vertices2 = vertices
                vertices1_coord, vertices2_coord = old_vor.vertices[[vertices1, vertices2]]
                polygon.append(vertices1_coord)
                polygon.append(vertices2_coord)
                # print(vertices1_coord,vertices1_coord)
            else:
                vertices1 = vertices[vertices != -1]
                vertices_1_coord = old_vor.vertices[vertices1]
                polygon.append(vertices_1_coord)
                mid_x, mid_y = 0.5 * (x + inter_p_coord[0]), 0.5 * (y + inter_p_coord[1])

                vertices_2_coord = modified_vertices_on_convex(dt, vertices_1_coord[0], vertices_1_coord[1], mid_x,
                                                               mid_y)
                polygon.append(vertices_2_coord)
                # print("vert",vertices1)
            polygons.append(polygon)
            # print("polygon",polygon)

            # print(new_vertices[vertices1],new_vertices[vertices2])
            # polygon.extend([new_vertices[vertices1],new_vertices[vertices2]])
    total_area = 0
    sum = 0
    # print(z_values,"asdfasd",polygons)
    for z, poly in zip(z_values, polygons):
        points = [[i[0], i[1], z] for i in poly]
        # print(points)
        area = polygon_area(points)
        total_area += area
        sum += z * area
        # print(area)
    Z = sum / total_area
    # print("Z value is ",Z)
    # start = time.time()
    # print("second roudn",start-end)
    z_hugo = dt.interpolate_nni(x, y)
    # print("z_hugo and z_max is",Z,z_hugo,x,y)
    if abs(Z - z_hugo) > 0.001:
        print("an erro is ", abs(Z - z_hugo))
    #print("z and z_hugo ", abs(Z - z_hugo))
    return Z

def find_surroundingPoint(dt,x,y):
    dt.insert_one_pt(x,y,1)
    vi = dt.closest_point(x, y)
    #print("triangular is",vi)
    tris =dt.adjacent_vertices_to_vertex(vi).tolist()
    #print(tris)

    #if 0 in tris:
    #   tris.remove(0)
    if len(tris) == 3:
        tris.append(0)
    pts = dt.points[tris]
    dt.remove(vi)
    return pts
    #tris.remove(vi)
