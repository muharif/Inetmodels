#!/usr/bin/env python
# coding: utf-8

import http.client
import json
import pandas as pd


class Inetmodels:

    def __init__(self, verbose = True):
        self.verbose = verbose
        self.addr = "inetmodels.com"
        self.header = {
            'content-type': "application/json",
#             'cache-control': "no-cache",
        }
        self.NT = {
            'GCN': "Gene Co-Expression Network",
            'MON': "Multi-Omics Network"
        }
        self.NetworkDict = self.__getNetworkDict()
        self.MON = dict(enumerate(self.Categories['MON'].keys()))
        self.GCN = dict(enumerate(self.Categories['GCN'].keys()))
        self.__printCategoryTypes()

        
    def printCategories(self, networkType = '', categoryType = 0):
        
        if networkType == 'MON':
            if ((categoryType >= 0) and (categoryType <= len(self.MON.keys()))):
                self.__printMONCategories(networkType, categoryType)
            else:
                raise ValueError('for networkType %s, categoryType must be between 0-%d' % (networkType,len(self.MON.keys())-1))

        elif networkType == 'GCN':
            if ((categoryType >= 0) and (categoryType <= len(self.GCN.keys()))):
                self.__printGCNCategories(networkType, categoryType)
            else:
                raise ValueError('for networkType %s, categoryType must be between 0-%d' % (networkType,len(self.GCN.keys())-1))        
        else:
            raise ValueError('networkType has to be MON or GCN, for Multi-Omics Network and Gene Co-Expression Network respectively')
    
    def query(self, 
              networkType = "", # MON or GCN, default = ""
              categoryType = 0,  # integer, default = 0
              categoryName = "",  # string, default = ""
              search = [], #list of strings, default = []
              pruning = 0, # #float, default = 2.5, means FDR < 10E-2.5
              nodeLimit = 10, #integer, default = 10 
              firstNeighbour = True, #boolean, default = True 
              gene = True, #boolean, default = True 
              metabolite = True, #boolean, default = True 
              protein = True, #boolean, default = True 
              clinical = True, #boolean, default = True 
              gut_microbiome = True, #boolean, default = True 
              oral_microbiome = True #boolean, default = True 
             ):
        
        #Sanity Check
        if networkType == 'MON':
            if ((categoryType >= 0) and (categoryType < len(self.MON.keys()))):
                if categoryName in self.Categories[networkType][self.MON[categoryType]].keys():
                    pass
                else:
                    raise ValueError("Invalid categoryName, check valid categoryName with 'printCategories('MON',%d)'" % categoryType)
            else:
                raise ValueError('for networkType %s, categoryType must be between 0-%d' % (networkType,len(self.MON.keys())-1))
        elif networkType == 'GCN':
            if ((categoryType >= 0) and (categoryType < len(self.GCN.keys()))):
                if categoryName in self.Categories[networkType][self.GCN[categoryType]].keys():
                    pass
                else:
                    raise ValueError("Invalid categoryName, check valid categoryName with 'printCategories('GCN',%d)'" % categoryType)
            else:
                raise ValueError('for networkType %s, categoryType must be between 0-%d' % (networkType,len(self.GCN.keys())-1))        
        else:
            raise ValueError('networkType has to be MON or GCN, for Multi-Omics Network and Gene Co-Expression Network respectively')
        
        
        #Check input type
        if type(networkType) != str:
            raise TypeError('networkType has to be string')
        if (type(categoryType) != int)  | (nodeLimit < 1):
            raise TypeError('categoryType has to be integer > 0')
        if type(categoryName) != str:
            raise TypeError('categoryName has to be string')
        if (type(search) != list) | (len(search) < 1):
            raise TypeError('search has to be list of analytes')
        if type(float(pruning)) != float:
            raise TypeError('pruning has to be float/decimal')
        if (type(nodeLimit) != int) | (nodeLimit < 1):
            raise TypeError('nodeLimit has to be integer > 0')
        if type(firstNeighbour) != bool:
            raise TypeError('firstNeighbour has to be True or False (boolean)')
        if type(gene) != bool:
            raise TypeError('gene has to be True or False (boolean)')
        if type(metabolite) != bool:
            raise TypeError('metabolite has to be True or False (boolean)')
        if type(protein) != bool:
            raise TypeError('protein has to be True or False (boolean)')       
        if type(clinical) != bool:
            raise TypeError('metabolite has to be True or False (boolean)')
        if type(oral_microbiome) != bool:
            raise TypeError('oral_microbiome has to be True or False (boolean)')   
        if type(gut_microbiome) != bool:
            raise TypeError('gut_microbiome has to be True or False (boolean)')  
            
        analytes = ',\n\t'.join(["\"%s\"" % i for i in search])+'\n'
        analyteTypes_lst = []
        if metabolite:
            analyteTypes_lst.append("METABOLITE")
        if protein:
            analyteTypes_lst.append("PROTEIN")
        if clinical:
            analyteTypes_lst.append("CLINICAL")
        if oral_microbiome:
            analyteTypes_lst.append("ORAL MICROBIOME")
        if gut_microbiome:
            analyteTypes_lst.append("GUT MICROBIOME")
        if gene:
            analyteTypes_lst.append("GENE")
        analyteTypes = ',\n\t'.join(["\"%s\"" % i for i in analyteTypes_lst])+'\n'
        if firstNeighbour:
            firstNeighbour = "true"
        else:
            firstNeighbour = "false"
        
        if networkType == 'MON':
            self.__queryMON(networkType,categoryType,categoryName,analytes,analyteTypes,pruning,nodeLimit,firstNeighbour)
        if networkType == 'GCN':
            self.__queryGCN(networkType, categoryType,categoryName,analytes,analyteTypes,pruning,nodeLimit,firstNeighbour)
        
        
    
    def __getNetworkDict(self):
        conn = http.client.HTTPSConnection(self.addr)
        conn.request("GET", "/api/data-types", "")
        self.Categories = json.loads(conn.getresponse().read())['categoryDataTypeMap']
        self.Categories['GCN'] = self.Categories.pop('Gene Co-Expression Network')
        self.Categories['MON'] = self.Categories.pop('Multi-Omics Network')
        conn.close()

    def __sendQuery(self):
        conn = http.client.HTTPSConnection(self.addr)
        conn.request("POST", "/api/query", self.__JSONquery, self.header)
        res = conn.getresponse()
        val = res.read()
        data = json.loads(val)
        conn.close()
        if len(data['links']) == 0:
            raise ValueError('No network has been found that fulfills the selected filters, please adjust it!')  
        edges = pd.DataFrame(data['links'])[['source','target','correlation','pvalue','padj']]
        nodes = pd.DataFrame(data['nodes']).set_index('id')
        id2index = nodes['index'].to_dict()
        edges['source'] = edges['source'].replace(id2index)
        edges['target'] = edges['target'].replace(id2index)
        edges.columns = ['Source','Target', 'Weight', 'P-Value', 'FDR']
        nodes_col = ['index' ,'symbol', 'info1', 'info2', 'info3', 'location']
        diff = set(nodes_col).difference(nodes.columns)
        for i in diff:
            nodes[i] = float('NaN')
        nodes = nodes[['index' ,'symbol', 'info1', 'info2', 'info3', 'location']]
        nodes.columns = ['Node','Symbol','Info1','Info2','Info3','AnalyteType']
        self.edges = edges
        self.nodes = nodes.set_index('Node')
        
    def __printCategoryTypes(self):
        print('networkType: MON --> Multi-Omics Network')
        for i in self.MON.keys():
            print('\tcategoryType %d: %s' % (i,self.MON[i]))
        print('networkType: GCN --> Gene Co-Expression Network')
        for i in self.GCN.keys():
            print('\tcategoryType %d: %s' % (i,self.GCN[i]))
        
    
    def __printMONCategories(self, networkType, categoryType = 1):
        print(networkType + ' (Multi-Omics Network)')
        print('\t%d: %s' % (categoryType,self.MON[categoryType]))
        for i in sorted(self.Categories[networkType][self.MON[categoryType]].keys()):
            print('\t\t' + i)
    
    def __printGCNCategories(self, networkType, categoryType = 1):
        print(networkType + ' (Gene Co-Expression Network)')
        print('\t%d: %s' % (categoryType,self.GCN[categoryType]))
        for i in sorted(self.Categories[networkType][self.GCN[categoryType]].keys()):
            print('\t\t' + i)
        
    
    def __payload(self, networkType,categoryType,categoryName,analytes,analyteTypes,pruning,nodeLimit,firstNeighbour):
        networkType = self.NT[networkType]
        self.__JSONquery = "{\n    \"networkType\": \"%s\",\n    \"categoryType\": \"%s\",\n    \"categoryName\": \"%s\",\n    \"analytes\" : [\n        %s    ],\n    \"analyteTypes\" : [\n        %s    ],\n    \"pruning\" : %d,\n    \"nodeLimit\" : %d,\n    \"firstNeighbour\" : %s\n}" % (networkType,categoryType,categoryName,analytes,analyteTypes,pruning,nodeLimit,firstNeighbour)
        if self.verbose:
            print('Your Request:')
            print(self.__JSONquery)
        
    def __queryGCN(self, networkType, categoryType,categoryName,analytes,analyteTypes,pruning,nodeLimit,firstNeighbour):
        categoryType = self.GCN[categoryType]
        self.__payload(networkType,categoryType,categoryName,analytes,analyteTypes,pruning,nodeLimit,firstNeighbour)
        self.__sendQuery()
        
    def __queryMON(self, networkType, categoryType,categoryName,analytes,analyteTypes,pruning,nodeLimit,firstNeighbour):
        categoryType = self.MON[categoryType]
        self.__payload(networkType,categoryType,categoryName,analytes,analyteTypes,pruning,nodeLimit,firstNeighbour)
        self.__sendQuery()