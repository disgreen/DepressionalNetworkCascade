# -*- coding: utf-8 -*-
"""
DepressionalNetwork (v11142022)
This is a module for creating a network of flow connected surface depressions, and simulating flow cascades between
them.
This version contains two classes: Depression and DepressionalNetwork
"""

import pandas as pd

class Depression(object):
    '''This class describes the depression object as described in Green and Crumpton (2022)
        self.maximum_volume = Vmax
        self.maximum_area = Amax
        self.maximum_depth = Hmax
        self.local_contributing_area = Acat + Amax
        self.local_catchment_area = Acat + Amax - Amax
        self.depression_type = Depression Type
    '''
    def __init__(self, identifier, maximum_volume, maximum_area, maximum_depth, local_contributing_area, upstream_depressions, downstream_depressions, depression_type = None, curve_number = 100, initial_volume = 0.0):
        self.identifier = identifier
        self.maximum_volume = maximum_volume
        self.maximum_area = maximum_area
        self.maximum_depth = maximum_depth
        self.local_contributing_area = local_contributing_area
        self.local_catchment_area = self.local_contributing_area - self.maximum_area 
        self.upstream_depressions = upstream_depressions
        self.downstream_depressions = downstream_depressions
        self.depression_type = depression_type
        self.local_curvenumber = curve_number
        self.local_precipitation = 0.0
        self.local_runoff = 0.0
        self.export = 0.0
        self.initial_volume = initial_volume
        self.storage_volume = 0.0
        self._upstream_export = 0.0
        self._upslope_maximum_storage_volume = 0.0
        self._upslope_contributing_area = 0.0
        self._upslope_storage_volume = 0.0
        self._runoff_contributing_area = 0.0
        self._runoff_contributing_area2 = 0.0
        self._upslope_specific_storage = 0.0
        self._number_upslope_contributing_depressions = 0
        self._strahler_number = 1
        self._shreve_number = 1
        self._runoff_network = set()
        self._upslope_network = set()
        self._downslope_network = set()
        self._outlet_depression = set()
    
    @property
    def runoff_network(self):
        '''
        This method gives the entire set of depressions that are upstream of self that have
        contributed runoff to self
        '''
        network = set()
        
        for depression in self.upstream_depressions:
            if depression.export > 0.0:
                if not depression in network:
                    network.add(depression)
                    network.update(depression.runoff_network)
        self._runoff_network = network
        return network
      
    @property
    def upslope_network(self):
        '''
        This method gives the entire set of depressions that are upstream of self
        '''
        network = set()
        self._upslope_network.clear()
            
        for depression in self.upstream_depressions:
            network.add(depression)
            network.update(depression.upslope_network)
        self._upslope_network = network
        return network
    
    @property
    def downslope_network(self):
        '''
        This method gives the entire set of depressions that are downslope of self
        '''
        network = set()
        self._downslope_network.clear()
        
        for depression in self.downstream_depressions:
            network.add(depression)
            network.update(depression.downslope_network)
        self._downslope_network = network
        return network
    
    @property
    def outlet_depression(self):
        '''
        This property identifies and retreives the outlet depression of the network the depression belongs to.
        '''
        network = set()
        self._outlet_depression.clear()
        if len(self.downslope_network) > 0:
            for depression in self.downslope_network:
                if depression.depression_type == 'Outlet':
                    network.add(depression)
        else:
            network.add(self)
        self._outlet_depression = network
        return network
    
    @property
    def upslope_specific_storage(self):
        '''
        This property calculates and retreives the total depressional specific storage upslope of the depression (not just adjacent neighboring features)
        '''
        sd = 0.0
        network = self.upslope_network
        for depression in network:
            sd += depression.maximum_volume/depression.local_contributing_area
        self._upslope_specific_storage = sd
        return sd
    
    @property
    def strahler_number(self):
        '''
        This property calculates the Strahler order of the depression
        '''
        order = 0
        nus = len(self.upstream_depressions)
        if nus == 0:
            order = 1
        else:
            sns = set()
            for depression in self.upstream_depressions:
                sns.add(depression.strahler_number)
    
            max_order = max(sns)
                   
            if len(sns) == 1:
                if nus > 1:
                    order = max_order + 1
                else:
                    order = max_order
            else:
                order = max_order
        
        self._strahler_number = order
        return self._strahler_number
    
    @property
    def shreve_number(self):
        '''
        This property calculates the Shreve order of the depression
        '''
        order = 0
        nus = len(self.upstream_depressions)
        if nus == 0:
            order = 1
        else:
            sns = []
            for depression in self.upstream_depressions:
                sns.append(depression.shreve_number)
    
            order = sum(sns)
            
        self._shreve_number = order
        return self._shreve_number
    
    @property
    def runoff_contributing_area(self):
        '''
        This property claculates the contributing area of all depressions contributing runoff to the depression, including but not limited to neighboring depressions
        '''
        area = 0.0
        network = self.runoff_network
        for depression in network:
            area += depression.local_contributing_area
        self._runoff_contributing_area = area
        return area
    
    @property
    def upslope_contributing_area(self):
        '''
        This property calculates the total contributing areas of all upslope depressions, including but not limited to neighboring features
        '''
        area = 0.0
        network = self.upslope_network
        for depression in network:
            area += depression.local_contributing_area
        self._upslope_contributing_area = area
        return area
    
    @property
    def upslope_storage_volume(self):
        '''
        This property claculates the total stored water volume of all depressions contributing runoff to the depression, including but not limited to neighboring depressions
        '''
        volume = 0.0
        network = self.upslope_network
        for depression in network:
            volume += depression.storage_volume
        self._upslope_storage_volume = volume
        return volume
    
    @property
    def upslope_maximum_storage_volume(self):
        '''
        This property claculates the total maximum water storage of all depressions contributing runoff to the depression, including but not limited to neighboring depressions
        '''
        volume = 0.0
        network = self.upslope_network
        for depression in network:
            volume += depression.maximum_volume
        self._upslope_maximum_storage_volume = volume
        return volume

    def calculate_local_runoff(self,precipitation):
        '''
        This method calculates local runoff using the NRCS curvenumber method
        Inputs: precipitation depth
        Outputs: runoff volume
        The method uses the average curve number assigned to the depression. 
        '''
        s = float(25.4/float(self.local_curvenumber) - 0.254)
        if precipitation <= 0.2*s:
            q = 0.0
        else:
            q = ((precipitation-0.2*s)**2)/(precipitation+0.8*s)
        local_runoff = q*self.local_contributing_area
        return local_runoff

    def accumulate_runoff_uniform(self,precipitation):
        '''
        This method routes runoff through the subnetwork upslope of the depression for a given spatially uniform precipitation depth
        '''
        upstream_export = 0.0
        self.export = 0.0
        self._runoff_contributing_area = 0.0
        self.storage_volume = self.initial_volume
        self.local_runoff = self.calculate_local_runoff(precipitation)
        self._runoff_network.clear()
        n = len(self.downstream_depressions)
        nds = (1 if n == 0 else n)
        
        if len(self.upstream_depressions) == 0:
            self.export = max([self.local_runoff - (self.maximum_volume-self.initial_volume),0.0])/float(nds)
            upstream_export = self.export
            self._upstream_export = self.local_runoff
            
            if self.export == 0.0:
                self.storage_volume += (self.local_runoff + self.initial_volume)
            else:
                self.storage_volume = self.maximum_volume 

        else:
            for depression in self.upstream_depressions:
                upstream_export += depression.accumulate_runoff_uniform(precipitation)
            
            self._upstream_export = upstream_export
            self.export = max([self.local_runoff + upstream_export - (self.maximum_volume-self.initial_volume),0.0])/float(nds)
            
            if self.export > 0.0:
                self.storage_volume = self.maximum_volume
                upstream_export += self.export
            else:
                self.storage_volume += (self.local_runoff + upstream_export + self.initial_volume)
        
        return self.export

    def __str__(self):
        obj_str = 'Identifier: {}\nDepression Type: {}\nMaximum Volume: {}\nMaximum Area: {}\nMaximum Depth: {}\nLocal Contributing Area: {}\nUpstream Depressions: {}\nDownstream Depressions: {}\nNRCS Curve Number: {}\nLocal Runoff: {}\nLocal Export: {}\nStorage Volume: {}\nTotal Upslope Contributing Area (Excludes self): {}\nUpslope Runoff Contributing Area (Excludes self): {}\nUpstream Contributing Depressions: {}'
        return obj_str.format(self.identifier,self.depression_type,self.maximum_volume,self.maximum_area,self.maximum_depth,self.local_contributing_area,
                               [d.identifier for d in self.upstream_depressions],[d.identifier for d in self.downstream_depressions],self.local_curvenumber,self.local_runoff,
                               self.export,self.storage_volume,self.upslope_contributing_area, self.runoff_contributing_area, [d.identifier for d in self.runoff_network])

class DepressionalNetwork:
  
    def __init__(self,name):
        self.name = name
        self.network = {}
        self._network_size = None
        self._network_list = []
        self._total_storage_volume = None
        self._total_contributing_area = None
        self._total_runoff_stored = None
        self._number_filled_depressions = None
        self._number_exporting_depressions = None
        self._outlet_depressions = []
        self._interior_depressions = []
        self._headwater_depressions = []
        self._disjunct_depressions = []
        self._runoff_contributing_depressions = set()
    
    @property
    def network_list(self):
        '''
        returns a list of depressions in the network
        '''
        if len(self._network_list) == 0:
            self._network_list = list(self.network.values())
        return self._network_list
    
    @property
    def number_of_depressions(self):
        """
        number_of_depressions (property) return the number of depressions in the network
        """
        self._network_size = len(self.network)
        return self._network_size
    
    @property
    def total_storage_volume(self):
        """
        total_storage_volume (property) return the sum of the maximum storage volumes of each depression
        """
        _network = self.network_list
        volume = 0.0
        for depression in _network:
            volume += depression.maximum_volume
        self._total_storage_volume = volume
        return self._total_storage_volume
    
    @property
    def total_runoff_stored(self):
        """
        total_runoff_stored (property) return the sum of the volume of runoff stored in each depression
        """
        _network = self.network_list
        storage = 0.0
        for depression in _network:
            storage += depression.storage_volume
        self._total_runoff_stored = storage
        return self._total_runoff_stored
    
    @property
    def total_contributing_area(self):
        """
        total_contributing_area (property) return the sum of all contributing area in the network
        """
        _network = self.network_list
        contributing_area = 0.0
        for depression in _network:
            contributing_area += depression.local_contributing_area
        self._total_contributing_area = contributing_area
        return self._total_contributing_area
    
    @property
    def number_of_depressions_filled(self):
        """
        number_of_depressions_filled (property) returns the number of depressions with volume >= maximum volume
        """
        _network = self.network_list
        number_filled = 0
        for depression in _network:
            if depression.storage_volume >= depression.maximum_volume:
                number_filled += 1
        self._number_of_depressions_filled = number_filled
        return self._number_of_depressions_filled
    
    @property
    def number_of_exporting_depressions(self):
        """
        number_of_exporting_depressions (property) returns the number of depressions with export > 0.0
        """
        _network = self.network_list
        number_exporting = 0
        for depression in _network:
            if depression.export > 0:
                number_exporting += 1
        self._number_of_exporting_depressions = number_exporting
        return self._number_of_exporting_depressions
    
    @property
    def outlet_depressions(self):
        """
        outlet_depressions (property) returns the set of depressions with type "Outlet"
        """
        if len(self._outlet_depressions) == 0:
            _network = self.network_list
            outlets = []
            for depression in _network:
                if depression.depression_type == 'Outlet':
                    outlets.append(depression)
            self._outlet_depressions = outlets
        return self._outlet_depressions
    
    @property
    def headwater_depressions(self):
        """
        headwater_depressions (property) returns the set of depressions with type "Headwater"
        """
        if len(self._headwater_depressions) == 0:
            _network = self.network_list
            headwaters = []
            for depression in _network:
                if depression.depression_type == 'Headwater':
                    headwaters.append(depression)
            self._headwater_depressions = headwaters
        return self._headwater_depressions
    
    @property
    def interior_depressions(self):
        """
        interior_depressions (property) returns the set of depressions with type "Interior"
        """
        if len(self._interior_depressions) == 0:
            _network = self.network_list
            interiors = []
            for depression in _network:
                if depression.depression_type == 'Interior':
                    interiors.append(depression)
            self._interior_depressions = interiors
        return self._interior_depressions
    
    @property
    def disjunct_depressions(self):
        """
        isolated_depressions (property) returns the set of depressions with type "Isolated" (refered to as disjunct in the manuscript).
        """
        if len(self._disjunct_depressions) == 0:
            _network = self.network_list
            disjuncts = []
            for depression in _network:
                if depression.depression_type == 'Disjunct':
                    disjuncts.append(depression)
            self._disjunct_depressions = disjuncts
        return self._disjunct_depressions
    
    def route_runoff_uniform_precipitation(self, precipitation):
        """
        route_runoff calls the accumulate_runoff (for outlet depressions) and calculate_local_runoff methods(for isolated depressions)
        inputs: self (depressional network), precipitation amount
        """
        Qr = 0.0
        Qt = 0.0
        Ar = 0.0        
        Dr = set()
        
        #Create collections of outlet and isolated depressions
        outlets = self.outlet_depressions
        disjuncts = self.disjunct_depressions
        
        #Empty the upstream_contributing_depressions set
        for depression in self.network_list:
            depression._runoff_network.clear()
            depression._runoff_contributing_area = 0.0
            Qt += depression.calculate_local_runoff(precipitation)
        
        #Loop over isolated and outlet depressions in the network and route runoff 
        for outlet in outlets:
            Qr += outlet.accumulate_runoff_uniform(precipitation) 
            if outlet.export > 0.0:
                Dr.update(outlet.runoff_network)
                Dr.add(outlet)
        
        for disjunct in disjuncts:
            Qr += disjunct.accumulate_runoff_uniform(precipitation) 
            if disjunct.export > 0.0:
                Dr.add(isolated)
        
        for dep in Dr:
            Ar += dep.local_contributing_area
        
        Nr = len(Dr)
        Nf = self.number_of_depressions_filled
        Vs = self.total_runoff_stored
        
        return Qr,Ar,Nr,Nf,Vs,Qt,Dr
            
    def construct_network(self, network_table):
        """
        construct_network (method) build the network set from the supplied pandas dataframe network table
        the network table should be constructed with the ConstructUpstreamSets and ConstructDownstreamSets functions
        prior to developing the cascade network. This version of the construct_network assumes that the US and DS sets are 
        given as python lists of integers.
        """
        if len(network_table)>0:
            for row in network_table.itertuples():
                identifier = row.DepressionID
                maximum_volume = row.MaximumVolume
                maximum_area = row.MaximumArea
                maximum_depth = row.MaximumDepth
                contributing_area = row.MinimumContributingArea
                depression_type = row.DepressionType
                curve_number = row.CurveNumber
                depression = Depression(identifier,maximum_volume,maximum_area,maximum_depth,contributing_area,[],[],depression_type,curve_number,0.0)
                self.network[identifier] = depression
                
            for identifier in self.network:
                upstream_depressions = network_table[network_table['DepressionID'] == identifier]['UpstreamDepressions'].values[0]
                downstream_depressions = network_table.loc[network_table['DepressionID'] == identifier]['DownstreamDepressions'].values[0]
                depression = self.network[identifier]
                
                if len(upstream_depressions) > 0:
                    for depid in upstream_depressions:
                        usdep = self.network[depid]
                        depression.upstream_depressions.append(usdep)
                
                if len(downstream_depressions) > 0:
                    for depid in downstream_depressions:
                        dsdep = self.network[depid]
                        depression.downstream_depressions.append(dsdep)
        else:
            self.network = {} 
        return self.network

#####Functions for constructing depressional networks within a Pandas dataframe, sourced from a text file
def ConstructDownstreamSets(deptable,jointable):
    '''This function constructs the set of downstream depressions for each depression in the network'''
    deptable['DownstreamDepressions'] = pd.Series(dtype='object')
    
    for row in deptable.itertuples():
        cdep = row.DepressionID
        indexset = jointable[jointable['DepressionID'] == cdep]
        
        depset = list()
        
        for subrow in indexset.itertuples():
            if subrow.CatchmentID != cdep:
                depset.append(int(subrow.CatchmentID))
        
        deptable.at[row.Index, 'DownstreamDepressions'] = depset

def ConstructUpstreamSets(deptable,jointable):
    '''This function constructs the set of upstream depressions for each depression in the network'''
    deptable['UpstreamDepressions'] = pd.Series(dtype='object')
    
    for row in deptable.itertuples():
        cdep = row.DepressionID
        
        indexset = jointable[jointable['CatchmentID'] == cdep]

        depset = list()
        
        for subrow in indexset.itertuples():
            if subrow.DepressionID != cdep:
                depset.append(int(subrow.DepressionID))

        deptable.at[row.Index, 'UpstreamDepressions'] = depset

def UpdateTypesAndFlags(deptable):
    '''This function classifies depressions according to their network type. 
    In the manuscript, the term disjunct' is used. 
    This is synonymous with 'isolated' in the code'''
    deptable['DuplicateFlag'] = pd.Series(dtype='str')
    deptable['DepressionType'] = pd.Series(dtype='str')
    
    for row in deptable.itertuples():
        usset = row.UpstreamDepressions
        dsset = row.DownstreamDepressions
        if (len(usset) == 0 and len(dsset) == 0):
            deptable.at[row.Index,'DepressionType'] = 'Disjunct'
            continue
        elif (len(usset) == 0 and len(dsset) > 0):
            deptable.at[row.Index,'DepressionType'] = 'Headwater'
            continue
        elif (len(usset) > 0 and len(dsset) == 0):
            deptable.at[row.Index,'DepressionType'] = 'Outlet'
            continue
        else:
            compset = set(usset) & set(dsset)
            if len(compset) > 0:
                deptable.at[row.Index,'DuplicateFlag'] = 'Duplicate'
            deptable.at[row.Index,'DepressionType'] = 'Interior'
            continue

def CheckNetworkStructure(deptable):
    '''This function checks the network structure for circular loops. 
    All depressions should route to the outlet feature eventually.'''
    deptable['StructureFlag'] = pd.Series(dtype='str')
    flag1 = ''
    flag2 = ''
    flag3 = ''
    flag4 = ''
    
    for row in deptable.itertuples():
        dep = row.DepressionID
        uset = row.UpstreamDepressions
        dset = row.DownstreamDepressions
        
        for _dep in uset:
            record = deptable[deptable['DepressionID'] == _dep]
            _dset = record['DownstreamDepressions'].values[0]
            index = record.index[0]
            
            if not dep in _dset:
                deptable.at[row.Index,'StructureFlag'] = flag1 + '|' + str(_dep)
                deptable.at[index,'StructureFlag'] = flag2 + '|' + str(dep)
        
        for _dep in dset:
            record = deptable[deptable['DepressionID'] == _dep]
            print(_dep)
            _uset = record['UpstreamDepressions'].values[0]
            
            if not dep in _uset:
                deptable.at[row.Index,'StructureFlag'] = flag3 + '|' + str(_dep)
                deptable.at[index,'StructureFlag'] = flag4 + '|' + str(dep)

        
