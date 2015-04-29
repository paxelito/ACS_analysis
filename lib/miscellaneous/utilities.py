#!/usr/bin/env python
# -*- coding: latin-1 -*-

'''
	This is a collection of useful script to adopt all over the project

'''

import sys, os

def zeroBeforeStrNum(tmpl, tmpL):
	''' 
		Function to create string zero string vector before graph filename.
		According to the total number of reactions N zeros will be add before the instant reaction number 
		(e.g. reaction 130 of 10000 the string became '00130')

		:param tmpl: number of digits of the number
		:param tmpL: number of digits of the MAX number
	'''
	
	strZero = ''
	nZeros = len(str(tmpL)) - len(str(tmpl))
	if nZeros > 0:
		for i in range(0,nZeros): strZero = strZero + '0'
	return strZero