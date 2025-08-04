""" Project: PfemApplication
    Developer: LMonforte
    Maintainer: LM
"""
# Futures
from __future__ import print_function, absolute_import, division

# Built-in/Generic Imports
import os.path
import numpy as np

# Kratos Imports
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import KratosMultiphysics.PfemApplication as KratosPfem
import KratosMultiphysics.ContactMechanicsApplication as KratosContact


def Factory(custom_settings, Model):
    if(type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string")
    return SumForces_anchor_new(Model, custom_settings["Parameters"])


class SumForces_anchor_new(KratosMultiphysics.Process):
    #

    def __init__(self, Model, custom_settings):

        KratosMultiphysics.Process.__init__(self)

        # settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
        "anchor_depth"        : 3.1,
        "velocity"            : 0.1,
        "model_part_name"     : "Main_Domain",
        "search_distance"     : 0.01,
        "length_shaft"        : 1.2,
        "length_tip"          : 0.3
        }
        """)

        settings = custom_settings
        settings.ValidateAndAssignDefaults(default_settings)
        
        self.model_part_name = settings["model_part_name"].GetString() 
        self.model = Model

        self.anchordepth = settings["anchor_depth"].GetDouble()
        self.lengthshaft = settings["length_shaft"].GetDouble()
        self.lengthtip = settings["length_tip"].GetDouble()
        self.Vy = settings["velocity"].GetDouble()
        
        self.searchdistance = settings["search_distance"].GetDouble()
        

        problem_path = os.getcwd()
        self.figure_path = os.path.join(
            problem_path, "forces_new.csv")

    def ExecuteBeforeOutputStep(self):
    
        self.model_part = self.model[self.model_part_name]

        delta_time = self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

        for node in self.model_part.GetNodes(0):
            delta_disp = node.GetSolutionStepValue(
                KratosMultiphysics.DISPLACEMENT)
            delta_disp = delta_disp - \
                node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 1)
            for ii in range(0, 2):
                delta_disp[ii] = delta_disp[ii] / delta_time
            if (node.SolutionStepsDataHas(KratosMultiphysics.VELOCITY)):
                node.SetSolutionStepValue(
                    KratosMultiphysics.VELOCITY, delta_disp)



    def GetPorePressureAnchor(self, XSearchBottom,XSearchTop,YSearchBottom,YSearchTop,SearchDistance):

        nodes = self.model_part.GetNodes()
        for node in nodes:
            if (node.HasDofFor(KratosMultiphysics.WATER_PRESSURE)):
                variable = KratosMultiphysics.WATER_PRESSURE
            elif (node.HasDofFor(KratosMultiphysics.PRESSURE)):
                variable = KratosMultiphysics.PRESSURE
            else:
                return 0.0


        vec_nodes_okay = []
        vec_nodes_okay_ids = []
       # YSearchTop = y_upper_limit
       # YSearchBottom = y_lower_limit
       # XSearchBottom=x_lower_anchor
       # XSearchTop=x_corner
        
        
        for node in self.model_part.GetNodes():
            Force = node.GetSolutionStepValue(KratosMultiphysics.CONTACT_FORCE)
            condition1 = abs(Force[0]) + abs(Force[1]) > 1e-8
            Normal = node.GetSolutionStepValue(KratosMultiphysics.NORMAL)
            condition2 = abs(Normal[0]) + abs(Normal[1]) > 1e-8 #and node.X < 3*self.radius)
            condition3 = abs(node.GetSolutionStepValue(variable)) > 1e-8
            condition3 = node.Is(KratosMultiphysics.RIGID) == False

            tan_angle_anchor = 0
            threshold_pos_x = XSearchBottom+SearchDistance
            # tangent of the inclination of the beam / anchor / whatever
            if abs(XSearchTop-XSearchBottom)>1e-8:
                tan_angle_anchor = (YSearchTop-YSearchBottom)/abs(XSearchTop-XSearchBottom)
                threshold_pos_x += (YSearchTop-node.Y)/tan_angle_anchor


            # define the maximum position along X for the node to be okay, depending on its vertical position and the angle of the anchor


            # use these above as conditions
            condition_mp1 = node.Y >= YSearchBottom and node.Y<= YSearchTop
            condition_mp2 = node.X <= threshold_pos_x

            
            
            # append the nodes in the desired location to a list
            if ((condition1 or condition2) and condition3):
                if (condition_mp1 and condition_mp2):
                    vec_nodes_okay.append(node)
                    vec_nodes_okay_ids.append(node.Id)
     #               print('node_x ' +str(node.X)+', thrs_x ' +str(threshold_pos_x)+', node_y '+str(node.Y)+ ' YB ' +str(YSearchBottom))
    #    print('Length vec nodes: '+str(len(vec_nodes_okay)))
        # Now i want to kwnow if there is an element that has these nodes
        total_pressure = 0
        nodes_counted = 0
        for elem in self.model_part.GetElements(0):
            a = [1, 2, 3]
            conec = []
            for thisNodeGeom in elem.GetNodes():
                thisNode = thisNodeGeom.Id
                if (thisNode in vec_nodes_okay_ids):
                        # if the node is part of an element (i.e. it is okay) and it's part of the region of space that you want -> sum its pressure to the total one and divide by the total number of nodes
                        ii = vec_nodes_okay_ids.index(thisNodeGeom.Id)
                        total_pressure+=vec_nodes_okay[ii].GetSolutionStepValue(variable)
                        nodes_counted+=1

        return total_pressure/nodes_counted
        
        
        



    def ExecuteAfterOutputStep(self):
        self.model_part = self.model[self.model_part_name]

        for node in self.model_part.GetNodes(0):

            if (node.SolutionStepsDataHas(KratosContact.CONTACT_STRESS)):
                CF = node.GetSolutionStepValue(KratosContact.CONTACT_STRESS)
                CF = 0.0*CF
                node.SetSolutionStepValue(KratosContact.CONTACT_STRESS, CF)

            if (node.SolutionStepsDataHas(KratosContact.EFFECTIVE_CONTACT_FORCE)):
                CF = node.GetSolutionStepValue(
                    KratosContact.EFFECTIVE_CONTACT_FORCE)
                CF = 0.0*CF
                node.SetSolutionStepValue(
                    KratosContact.EFFECTIVE_CONTACT_FORCE, CF)

            if (node.SolutionStepsDataHas(KratosContact.EFFECTIVE_CONTACT_STRESS)):
                CF = node.GetSolutionStepValue(
                    KratosContact.EFFECTIVE_CONTACT_STRESS)
                CF = 0.0*CF
                node.SetSolutionStepValue(
                    KratosContact.EFFECTIVE_CONTACT_STRESS, CF)

            if (node.SolutionStepsDataHas(KratosMultiphysics.CONTACT_NORMAL)):
                CF = node.GetSolutionStepValue(
                    KratosMultiphysics.CONTACT_NORMAL)
                CF = 0.0*CF
                node.SetSolutionStepValue(
                    KratosMultiphysics.CONTACT_NORMAL, CF)

            if (node.SolutionStepsDataHas(KratosSolid.DISPLACEMENT_REACTION)):
                # print(node)
                for step in range(0, 2):
                    CF = node.GetSolutionStepValue(
                        KratosSolid.DISPLACEMENT_REACTION, step)
                    normal = node.GetSolutionStepValue(
                        KratosMultiphysics.NORMAL, step)
                    #print('NODE', node.Id, ' normal', normal)
                    #print('  step ', step, 'REACTION', CF)
                    CF = 0.0*CF
                    #print('  step ', step, 'REACTION', CF)
                    node.SetSolutionStepValue(
                        KratosSolid.DISPLACEMENT_REACTION, step, CF)
                # print('----------')

    #
    def ExecuteFinalizeSolutionStep(self):
        self.model_part = self.model[self.model_part_name]
        Q = self._GetResistance()

        Q2 = self._GetUpperForce()
        
        Q3 = self._GetResistance_x()
        
        Q4 = self._GetUpperForce_x()

        time = self._GetStepTime()
        
        funz1=self.Vy*(self._GetStepTime()-0.1)
        funz2=self.Vy/10*(self._GetStepTime()*np.pi/0.2-np.sin(self._GetStepTime()*np.pi/0.2))/np.pi
        
        corner_position = self.anchordepth + max(funz1,funz2)
        search_distance=self.searchdistance
        SearchDistance = search_distance
        y_upper_limit = corner_position + self.lengthshaft
        y_lower_limit = corner_position - self.lengthtip
        x_corner = 0.11
        x_lower_anchor = 0.125
        #print('CORNER:'+str(corner_position)+'\n')
        #print('SD:'+str(search_distance)+'\n')
              
        
        U_shaft = self.GetPorePressureAnchor(x_corner,x_corner,corner_position,y_upper_limit,SearchDistance)
        
        print(U_shaft)
        
        U_tip = self.GetPorePressureAnchor(x_lower_anchor,x_corner,y_lower_limit,corner_position,search_distance)
      
        print('U_tip:' +str(U_tip)+ ', U_shaft: ' +str(U_shaft))


        line_value = str(time) + " , " + str(corner_position) + " , " + str(Q)+ " , " + str(Q2) + " , " + str(Q3) + " , " + str(Q4) + " , " + str(U_shaft)+ " , " + str(U_tip)+  ",  \n"

        figure_file = open(self.figure_path, "a")
        figure_file.write(line_value)
        figure_file.close()



    # function to try to unfix the water pressure of nodes in that are in contact
    def TryToUnFixNodes(self):
        
        funz1=self.Vy*(self._GetStepTime()-0.1)
        funz2=self.Vy/10*(self._GetStepTime()*np.pi/0.2-np.sin(self._GetStepTime()*np.pi/0.2))/np.pi
        YLim = self.Y0 + max(funz1,funz2)

        if (self.dissipation_set):
            if (YLim < self.dissipation_depth):
                YLim = self.dissipation_depth

        for node in self.model_part.GetNodes(0):
            Force = node.GetSolutionStepValue(KratosMultiphysics.CONTACT_FORCE)
            Normal = node.GetSolutionStepValue(KratosMultiphysics.NORMAL)
            IsBoundary = False
            if (abs(Normal[0]) + abs(Normal[1]) > 1.0e-7):
                IsBoundary = True

            if (node.HasDofFor(KratosMultiphysics.WATER_PRESSURE)):
                if ((abs(Force[0]) + abs(Force[1]) > 1.0e-7)):
                    if (node.IsFixed(KratosMultiphysics.WATER_PRESSURE)):
                        node.Free(KratosMultiphysics.WATER_PRESSURE)
                        print(' UNFIXING ')
                elif (node.Y > YLim and node.X < 4.0*self.radius and IsBoundary):
                    if (node.IsFixed(KratosMultiphysics.WATER_PRESSURE) == False):
                        node.Fix(KratosMultiphysics.WATER_PRESSURE)
                        print(' UNFIXING :: FIXING ')

                elif (node.Y > YLim and node.X < 4.0*self.radius and (IsBoundary == False)):
                    if (node.IsFixed(KratosMultiphysics.WATER_PRESSURE)):
                        print(' UNFIXING :: INNER NODE :: DISASTER ')
                        node.Free(KratosMultiphysics.WATER_PRESSURE)

    # make the other stupid file just to check everything is sort of correct
    def MakeTheOtherFile(self):

        time = self._GetStepTime()

        problem_path = os.getcwd()
        self.other_file = os.path.join(problem_path, "other_file_data.csv")
        self.other_file = open(self.other_file, "a")

        time = str(time) + " \n "
        self.other_file.write(time)

        YLim = self.Y0 - self.Vy * self._GetStepTime()
        YLim = self._GetU2Position()

        for node in self.model_part.GetNodes(0):
            Force = node.GetSolutionStepValue(KratosMultiphysics.CONTACT_FORCE)
            if ((abs(Force[0]) + abs(Force[1]) > 1.0e-7) or (node.X0 < 1e-5 and node.Y > YLim-0.1)):
                x = node.X
                y = node.Y
                Force = node.GetSolutionStepValue(
                    KratosMultiphysics.CONTACT_FORCE)
                Stress = node.GetSolutionStepValue(
                    KratosContact.CONTACT_STRESS)
                EfForce = node.GetSolutionStepValue(
                    KratosContact.EFFECTIVE_CONTACT_FORCE)
                EfStress = node.GetSolutionStepValue(
                    KratosContact.EFFECTIVE_CONTACT_STRESS)
                WP = node.GetSolutionStepValue(
                    KratosMultiphysics.WATER_PRESSURE)

                line = str(x) + "  " + str(y) + "  " + \
                    str(YLim) + " " + str(WP) + " "
                line = line + self.AddVector(Force)
                line = line + self.AddVector(Stress)
                line = line + self.AddVector(EfForce)
                line = line + self.AddVector(EfStress)
                line = line + "\n"
                self.other_file.write(line)

        self.other_file.close()

    def AddVector(self, vector):
        line = str(vector[0]) + " " + str(vector[1]) + " "
        return line
    #
    def _GetStepTime(self):

        return self.model_part.ProcessInfo[KratosMultiphysics.TIME]

    #
    def _GetStepDeltaTime(self):

        return self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

    def _GetU2Position(self):
        avance = self.Vy*self._GetStepTime()
        if (self.dissipation_set):
            if (avance > self.dissipation_depth):
                avance = self.dissipation_depth
        YSearch = self.Y0 - avance
        return YSearch

    def _GetPorePressureU2(self):
        YSearch = self._GetU2Position()
        U22 = self._GetPorePressureShaft(YSearch)
        return U22

    def _GetPorePressureU3(self):
        YSearch = self._GetU2Position()
        YSearch = YSearch + 7.5*self.radius
        U33 = self._GetPorePressureShaft(YSearch)
        return U33

    def _GetPorePressureU1(self):
        YSearch = self._GetU2Position()
        YSearch = YSearch - self.radius / np.tan(30.0*3.14159/180.0)/2.0
        U11 = self._GetPorePressureShaft(YSearch)
        return U11

    def _GetResistance(self):
        result = 0.0
        for node in self.model_part.GetNodes(0):
                Force = node.GetSolutionStepValue(
                    KratosMultiphysics.CONTACT_FORCE)
                result = result + Force[1]

        return result
        
    def _GetResistance_x(self):
        result = 0.0
        
        for node in self.model_part.GetNodes(0):
                Force = node.GetSolutionStepValue(
                    KratosMultiphysics.CONTACT_FORCE)
                result = result + Force[0]

        return result

    def _GetUpperForce(self):
        result = 0.0

        funz1=self.Vy*(self._GetStepTime()-0.1)
        funz2=self.Vy/10*(self._GetStepTime()*np.pi/0.2-np.sin(self._GetStepTime()*np.pi/0.2))/np.pi
        
        YMin = self.anchordepth + max(funz1,funz2)

        for node in self.model_part.GetNodes(0):
            if ( node.Y >= YMin):
                Force = node.GetSolutionStepValue(
                    KratosMultiphysics.CONTACT_FORCE)
                result = result + Force[1]

        return result
        
    def _GetUpperForce_x(self):
        result = 0.0

        funz1=self.Vy*(self._GetStepTime()-0.1)
        funz2=self.Vy/10*(self._GetStepTime()*np.pi/0.2-np.sin(self._GetStepTime()*np.pi/0.2))/np.pi
        
        YMin = self.anchordepth + max(funz1,funz2)

        for node in self.model_part.GetNodes(0):
            if ( node.Y >= YMin):
                Force = node.GetSolutionStepValue(
                    KratosMultiphysics.CONTACT_FORCE)
                result = result + Force[0]

        return result
    #
    def _GetFriction(self):
        result = 0
        YMin = self._GetU2Position()
        YMax = YMin + 7.5*self.radius

        for node in self.model_part.GetNodes(0):
            if (node.Y >= YMin):
                if (node.Y <= YMax):
                    Force = node.GetSolutionStepValue(
                        KratosMultiphysics.CONTACT_FORCE)
                    result = result + Force[1]

        return result

    def _GetPorePressureShaft(self, YSearch):

        nodes = self.model_part.GetNodes()
        for node in nodes:
            if (node.HasDofFor(KratosMultiphysics.WATER_PRESSURE)):
                variable = KratosMultiphysics.WATER_PRESSURE
            elif (node.HasDofFor(KratosMultiphysics.PRESSURE)):
                variable = KratosMultiphysics.PRESSURE
            else:
                return 0.0

        YBestTop = 100000000
        YBestBottom = -100000000
        nBestTop = -100
        nBestBottom = -100

        for node in self.model_part.GetNodes(0):
            Force = node.GetSolutionStepValue(KratosMultiphysics.CONTACT_FORCE)
            condition1 = abs(Force[0]) + abs(Force[1]) > 1e-8
            Normal = node.GetSolutionStepValue(KratosMultiphysics.NORMAL)
            condition2 = (abs(Normal[0]) + abs(Normal[1])
                          > 1e-8 and node.X < 3*self.radius)
            condition3 = abs(node.GetSolutionStepValue(variable)) > 1e-8
            condition3 = node.Is(KratosMultiphysics.RIGID) == False
            if ((condition1 or condition2) and condition3):
                YThis = node.Y
                if ((YThis-YSearch) > 0):
                    if (abs(YThis-YSearch) < abs(YBestTop - YSearch)):
                        nBestTop = node.Id
                        YBestTop = YThis
                elif ((YThis-YSearch) <= 0):
                    if (abs(YThis-YSearch) < abs(YBestBottom - YSearch)):
                        nBestBottom = node.Id
                        YBestBottom = YThis

        if ((nBestTop < 1) or (nBestBottom < 1)):
            print(
                ' In the Usomething. NotFound Contacting nodes that are in the range of this kind of thing')
            return 0.0

        # Now i want to kwnow if there is an element that has these nodes
        ReallyFound = False
        for elem in self.model_part.GetElements(0):
            a = [1, 2, 3]
            conec = []
            found = 0
            for thisNodeGeom in elem.GetNodes():
                thisNode = thisNodeGeom.Id
                if (thisNode == nBestTop):
                    found = found + 1
                if (thisNode == nBestBottom):
                    found = found + 1
            if (found == 2):
                ReallyFound = True
                break

        if (ReallyFound == False):
            print(' In the U Something. The two nodes do not share an element ')
            return 0.0

        DeltaY = abs(YBestBottom - YBestTop)
        NTop = 1 - abs(YBestTop - YSearch) / DeltaY
        NBottom = 1 - NTop

        if (NTop > 1.0 or NTop < 0):
            print('ULTRA MEGA STUPID ERROR ')

        if (NBottom > 1.0 or NBottom < 0):
            print('ULTRA MEGA STUPID ERROR ')

        uBottom = nodes[nBestBottom].GetSolutionStepValue(variable)
        uTop = nodes[nBestTop].GetSolutionStepValue(variable)
        ThisV = NTop*uTop + NBottom*uBottom
        return ThisV
