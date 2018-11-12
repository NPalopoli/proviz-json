#!/usr/bin/python

# proviz-json - Make JSON files for Proviz from ProtMiscuity data
# Copyright (C) 2018 Nicolas Palopoli
#
# # This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# # This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Module:       proviz-json
Description:  Make JSON files for Proviz from Protmiscuity data
Version:      0.2
Python:       2.7.12
Last Edit:    2018-11-10
Functions:
	None
Class:
	protmiscuity = data from ProtMiscuity
Objects:
	(protmiscuity)
    protmiscuityObj = main data from ProtMiscuity
Methods:
	(protmiscuity)
	check_data_directory(),check_output_directory(),parseSDT(),
	parseProtMiscuity(),mergeProtMiscuity(),findRanges(),makeJSONPerSiteData(),
	makeJSONPerSite(),makeProtMiscuityJSON()
Commandline:
    ./proviz-json.py
Dependencies:
	General: os,json,collections,more_itertools
"""

################################################################################
## SECTION I: GENERAL SETUP & PROGRAM DETAILS                                 ##
################################################################################

# General imports
import os
import json
from collections import OrderedDict
from more_itertools import consecutive_groups

################################################################################
## END OF SECTION I                                                           ##
################################################################################

################################################################################
## SECTION II: CLASSES                                                        ##
################################################################################

##----------------------------------------------------------------------------##
# protmiscuity class - data from ProtMiscuity
##----------------------------------------------------------------------------##
class protmiscuity():

    ##------------------------------------------------------------------------##
	# __init__ method - Define initial state of object
	##------------------------------------------------------------------------##
	def __init__(self):
		'''
		Define initial state of protmiscuity object
		'''
		# Configuration options
		self.options = {}
		self.options["base_path"] = ""
		self.options["data_path"] = ""
		self.options["output_path"] = ""
		# Data
		self.protein_data = {}
		self.protein_json = {}
		self.protein_fields = ("can_act_site","prom_act_site")

	##------------------------------------------------------------------------##
	# check_data_directory function - Ensure data path exists
	##------------------------------------------------------------------------##
	def check_data_directory(self):
		'''
		Ensure path exists
		'''
		self.options["data_path"] = os.path.join(self.options["base_path"], "tablas")
		if not os.path.exists(self.options["data_path"]):
			os.mkdir(self.options["data_path"])

	##------------------------------------------------------------------------##
	# check_output_directory function - Ensure output path exists
	##------------------------------------------------------------------------##
	def check_output_directory(self):
		'''
		Ensure output path exists
		'''
		self.options["output_path"] = os.path.join(self.options["base_path"], "proviz-json", "json")
		if not os.path.exists(self.options["output_path"]):
			os.mkdir(self.options["output_path"])

    ##------------------------------------------------------------------------##
	# parseSDT function - Parse semicolon-delimited text file
	##------------------------------------------------------------------------##
	def parseSDT(self, tsv_path, sdt_key):
		'''
		Parse tab-delimited text file into dict
		'''
		# Create empty dicts
		sdt_data = {}
		headers = []
		# Iterate over file
		for line in open(tsv_path).read().strip().split("\n"):
			if len(line) > 0:  # avoid empty lines
				lineBits = line.strip().replace('"', '').split(";")  # Remove quotes and split fields
				if line[0] == "#":  # avoid comments
					pass
				else:
					if headers == []:  # populate header if empty
						headers = lineBits
					else:
						try:
							# Create new entry in dict if key not present
							if lineBits[headers.index(sdt_key)] not in sdt_data:
								sdt_data[lineBits[headers.index(sdt_key)]] = []
							# Populate temp dict
							tmp_dict = {}
							for i in range(0, len(headers)):  # parse data according to headers
								if len(lineBits) > i:  # iterate while fields are present in line
									tmp_dict[headers[i]] = lineBits[i]  # set value of present fields
								else:
									tmp_dict[headers[i]] = "?"  # set value of missing fields
							# Append temp dict entry to list in output dict
							sdt_data[lineBits[headers.index(sdt_key)]].append(tmp_dict)
						# Skip lines if needed
						except Exception, e:
							print "skipping line", lineBits
							print e
		# Return dict of tab-delimited text
		# (keys: sdt_key values; values: list (each line an element) of dicts (each field a key:value pair))
		return sdt_data

    ##------------------------------------------------------------------------##
	# parseProtMiscuity function - Parse ProtMiscuity tables
	##------------------------------------------------------------------------##
	def parseProtMiscuity(self,table_name,key_field):
		'''
		Parse raw ProtMiscuity tables
		'''
		# Define path to ProtMiscuity protein data file
		sdt_path = os.path.join(protmiscuityObj.options["data_path"],table_name)
		# Parse colon-delimited ProtMiscuity data file to dictionary
		pm_proteina = protmiscuityObj.parseSDT(sdt_path,key_field)
		return(pm_proteina)

	##------------------------------------------------------------------------##
	# mergeProtMiscuity function - Merge parsed ProtMiscuity tables into dict
	##------------------------------------------------------------------------##
	def mergeProtMiscuity(self,pm_dict):
		'''
		Merge parsed ProtMiscuity tables into dict
		'''
		# Define list of relevant fields
		field_proteina = ["secuencia","codigo_uni_prot"]
		field_sitios_activos = ["can_act_site","prom_act_site"]
		# Identify table by use of relevant fields as dict keys
		if all (field in pm_dict[next(iter(pm_dict))][0] for field in field_proteina):
			field_names = field_proteina
		elif all (field in pm_dict[next(iter(pm_dict))][0] for field in field_sitios_activos):
			field_names = field_sitios_activos
		# Iterate over protein IDs
		for protein_id in pm_dict:
			# Create key in protein data dict if needed
			if not protein_id in self.protein_data:
				self.protein_data[protein_id] = {}
				for field in field_proteina:
					self.protein_data[protein_id][field] = ''
				for field in field_sitios_activos:
					self.protein_data[protein_id][field] = ''
			# Iterate over relevant fields
			for field in field_names:
				# Copy values into protein data dict
				self.protein_data[protein_id][field] = pm_dict[protein_id][0][field]
		return()

	##------------------------------------------------------------------------##
	# findRanges function - Get ranges of consecutive numbers in list
	##------------------------------------------------------------------------##
	def findRanges(self,iterable):
		'''
		Get ranges of consecutive numbers in list
		'''
		# Sort numbers
		iterable.sort()
		# Iterate over lists of consecutive numbers
	   	for group in consecutive_groups(iterable):
			group = list(group)
			# Return ranges
			if len(group) == 1:
				yield([group[0]])
			else:
				yield([group[0],group[-1]])

	##------------------------------------------------------------------------##
	# makeJSONPerSiteData function - Make JSON structure for data field
	##------------------------------------------------------------------------##
	def makeJSONPerSiteData(self,site_type,site_pos,sequence):
		'''
		Make JSON structure for data field in types of sites
		'''
		# Remove non-determined positions from list
		site_pos_list = site_pos.split(',')
		if 'ND' in site_pos_list:
			site_pos_list.remove('ND')
		# Make list of positions as int
		if '' in site_pos_list:
			pos_ranges = ['']
			site_data_list = []
		else:
			# Make list of position ranges
			pos_list = map(int,site_pos_list)
			site_data_list = []
			pos_ranges = list(self.findRanges(pos_list))
			# Iterate over positions
			for pos_range in pos_ranges:
				# Make ordered dict for data field
				site_data_dict = OrderedDict()
				# Get start and end positions in range
				site_data_dict["start"] = str(pos_range[0])
				site_data_dict["end"] = str(pos_range[-1])
				# Copy sequence in range
				site_data_dict["sequence"] = sequence[int(site_data_dict["start"])-1:int(site_data_dict["end"])]
				# Define colour
				site_data_dict["colour"] = "#0087ff"
				# Get type of site
				if site_type == "can_act_site":
					site_data_dict["hover"] = "canonic site"
				elif site_type == "prom_act_site":
					site_data_dict["hover"] = "promiscuous site"
				# Add fields to list in data field
				site_data_list.append(site_data_dict)
		# Return list of fields in data field
		return(site_data_list)

	##------------------------------------------------------------------------##
	# makeJSONPerSite function - Make JSON structure for types of sites
	##------------------------------------------------------------------------##
	def makeJSONPerSite(self,site_type,site_pos,sequence):
		'''
		Make JSON structure for types of sites
		'''
		# Make ordered dict for site data
		site_data = OrderedDict()
		# Define common fields
		site_data["type"] = "peptides"
		site_data["position"] = "-1"
		# Define fields specific of type of site
		if site_type == "can_act_site":  # canonic site
			site_data["name"] = "canonic site"
			site_data["colour"] = "#0087ff"
			site_data["help"] = "Canonic site residues annotated in ProtMiscuity"
			site_data["text_colour"] = "#000"
			# Make JSON structure for data field
			site_data["data"] = self.makeJSONPerSiteData(site_type,site_pos,sequence)
		elif site_type == "prom_act_site":  # promiscuous site
			site_data["name"] = "promiscuous site"
			site_data["colour"] = "#FF0054"
			site_data["help"] = "Promiscuous site residues annotated in ProtMiscuity"
			site_data["text_colour"] = "#000"
			# Make JSON structure for data field
			site_data["data"] = self.makeJSONPerSiteData(site_type,site_pos,sequence)
		# Return dict with data for all JSON fields
		return(site_data)

	##------------------------------------------------------------------------##
	# makeProtMiscuityJSON function - Make JSON files from merged tables
	##------------------------------------------------------------------------##
	def makeProtMiscuityJSON(self,protein_id,pm_dict):
		'''
		Make JSON files from merged ProtMiscuity tables
		'''
		# Make a dict of site types
		pm_sites = {field: pm_dict[field] for field in self.protein_fields}
		# Iterate over types of sites
		for site_type in pm_sites:
			# Make JSON structure only for site types present in protein
			if pm_sites[site_type]:
				self.protein_json[protein_id][site_type] = self.makeJSONPerSite(site_type,pm_sites[site_type],pm_dict["secuencia"])
		return()

################################################################################
## END OF SECTION II                                                          ##
################################################################################

################################################################################
## SECTION III: MAIN PROGRAM                                                  ##
################################################################################

if __name__ == "__main__":
	# Create object of class protmiscuity
	protmiscuityObj = protmiscuity()
	# Set paths
	protmiscuityObj.options["base_path"] = "/home/npalopoli/ProtMiscuity"
	protmiscuityObj.check_data_directory()
	protmiscuityObj.check_output_directory()
	# Parse ProtMiscuity protein data file
	pm_proteina = protmiscuityObj.parseProtMiscuity("proteina.csv","id")
	protmiscuityObj.mergeProtMiscuity(pm_proteina)
	# Parse ProtMiscuity active site data file
	pm_sitios_activos = protmiscuityObj.parseProtMiscuity("sitio_activo.csv","proteina_id")
	protmiscuityObj.mergeProtMiscuity(pm_sitios_activos)
	# Iterate over protein IDs
	for protein_id in protmiscuityObj.protein_data:
		# Create ordered dict for JSON data
		protmiscuityObj.protein_json[protein_id] = OrderedDict()
		# Make JSON files from merged ProtMiscuity tables
		protmiscuityObj.makeProtMiscuityJSON(protein_id,protmiscuityObj.protein_data[protein_id])
		# Make path to output JSON file
		json_path = os.path.join(protmiscuityObj.options["output_path"],''.join([protmiscuityObj.protein_data[protein_id]["codigo_uni_prot"],".json"]))
		# Make list of data for site types present in protein
		protein_json_list = [protmiscuityObj.protein_json[protein_id][field] for field in protmiscuityObj.protein_fields if field in protmiscuityObj.protein_json[protein_id]]
		# Write JSON file only with data for site types present in protein
		if protein_json_list:
			with open(json_path, 'w') as outfile:
				json.dump(protein_json_list, outfile)
	
################################################################################
## END OF SECTION III                                                         ##
################################################################################