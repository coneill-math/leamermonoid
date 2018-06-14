from sage.combinat.posets.posets import FinitePoset



class LeamerMonoid:
	
	def __init__(self, gens, s, init_internals = True):
		self.s = s
		
		#self.gapname = 'x'.join([str(i) for i in gens]) + 'x' + str(s) + 'xGammaInternal'
		self.gamma = NumericalSemigroup(gens)
		self.frob = self.gamma.FrobeniusNumber()
		self.gens = self.gamma.gens
		
		self.x0 = 0
		self.xf = 0
		self.heights = {}
		self.irrs = {}
		self.__lambda = 0
		
		self.factorizations = {}
		self.catenaries = {}
		
		if init_internals:
			self.FindIrreducibles()
	
	def HeightOfRow(self,x):
		self.FindHeights()
		
		if x > self.frob:
			return oo
		elif x <= 0:
			return 0
		else:
			return self.heights[x]
	
	def FirstIrreducibleAt(self,x):
		self.FindIrreducibles()
		
		if x > self.frob + self.x0:
			return oo
		else:
			return self.irrs[x]
	
	def LastReducibleAt(self,x):
		self.FindIrreducibles()
		
		if self.HeightOfRow(x) < oo and self.FirstIrreducibleAt(x) == oo:
			return self.HeightOfRow(x)
		else:
			return self.FirstIrreducibleAt(x)-1
	
	def IsInGamma(self, x, n):
		self.FindHeights()
		
		return self.HeightOfRow(x) >= n
	
	def IsIrreducible(self, x, n):
		self.FindIrreducibles()
		
		return self.IsInGamma(x,n) and (n == 1 or self.FirstIrreducibleAt(x) <= n)
	
	def Factorizations(self,x,n):
		self.FindIrreducibles()
		
		if (x,n) in self.factorizations.keys():
			return self.factorizations[(x,n)]
		
		if self.IsIrreducible(x,n):
			self.factorizations[(x,n)] = [[(x,n)]]
			return [[(x,n)]]
		
		facts = []
		
		for i in (1 .. (x-self.x0)):
			for j in (1 .. min(self.HeightOfRow(i),n)):
				if not(self.IsIrreducible(i,j)) or not(self.IsInGamma(x-i,n-j)):
					continue
				
				subfacts = self.Factorizations(x-i,n-j)
				for fact in subfacts:
					if fact[-1][0] < i or (fact[-1][0] == i and fact[-1][1] <= j):
						facts.append(fact + [(i,j)])
		
		self.factorizations[(x,n)] = facts
		
		return facts
	
	def FactorizationsOfLength(self,x,n,l):
		return [fact for fact in self.Factorizations(x,n) if len(fact) == l]
	
	def LengthSet(self,x,n):
		return Set([len(fact) for fact in self.Factorizations(x,n)])
	
	def LengthList(self,x,n):
		return sorted(list(Set([len(fact) for fact in self.Factorizations(x,n)])))
	
	def DeltaSet(self,x,n):
		lengths = sorted([i for i in self.LengthSet(x,n)])
		
		delta = []
		for i in (0 .. len(lengths) - 2):
			delta.append(lengths[i+1] - lengths[i])
		
		return Set(delta)
	
	def DeltaSetBelow(self, xmax = 0, ymax = 0):
		if xmax == 0:
			xmax = self.xmax + 2*self.x0
		if ymax == 0:
			ymax = self.ymax
		
		deltaset = Set([])
		
		for x in (1 .. xmax):
			for n in (1 .. min(self.HeightOfRow(x),ymax)):
				deltaset = deltaset.union(self.DeltaSet(x,n))
		
		return deltaset
	
	def DeltaSetCutoff(self):
		xmax = self.xmax + 2*self.x0
		ymax = self.ymax
		
		maxDelta = max(self.DeltaSetBelow())
		
		deltaSet = Set([])
		
		x = self.x0
		n = 1
		retarray = {i:0 for i in (1 .. maxDelta)}
		
		while len(deltaSet) < maxDelta:
			if self.HeightOfRow(x) == 0:
				x = x + 1
				continue
			
			for n in (1 .. min(self.HeightOfRow(x),ymax)):
				deltaSet = deltaSet.union(self.DeltaSet(x,n))
				for d in deltaSet:
					if retarray[d] == 0:
						retarray[d] = (x,n)
				
				if len(deltaSet) == maxDelta:
					break
			
			x = x + 1
		
		return retarray
	'''
	def DeltaSetBelowLarge(self, xmax = 0, ymax = 0):
		if xmax == 0:
			xmax = self.xmax + 2*self.x0
		if ymax == 0:
			ymax = self.ymax
		
		retlist = []
		
		for x in (1 .. xmax):
			for n in (1 .. min(self.HeightOfRow(x),ymax)):
				deltaset = self.DeltaSet(x,n)
				if len(deltaset) > 2 or (len(deltaset) > 1 and not(1 in deltaset)):
					retlist.append(((x,n), deltaset))
		
		return retlist
	'''
	def DeltaSetsWith(self, delta, xmax = 0, ymax = 0):
		if xmax == 0:
			xmax = self.xmax + 2*self.x0
		if ymax == 0:
			ymax = self.ymax
		
		retlist = []
		
		for x in (1 .. xmax):
			for n in (1 .. min(self.HeightOfRow(x),ymax)):
				if delta in self.DeltaSet(x,n):
					retlist.append((x,n))
		
		return retlist
	
	def LengthSetTable(self, xmin, xmax, ymax):
		T = [[0] + [0 for j in (xmin..xmax)] for i in (1 .. ymax)]
		for i in (1 .. ymax):
			T[i-1][0] = i
			for j in (xmin..xmax):
				T[i-1][j-xmin+1] = self.LengthSet(j,i)
		T = [[''] + [j for j in (xmin..xmax)]] + T
		return T
	
	def HTMLLengthSetTable(self, xmin, xmax, ymax):
		return table(self.LengthSetTable(xmin,xmax,ymax))
	
	def DeltaSetTable(self, xmin, xmax, ymax):
		T = [[0] + [0 for j in (xmin..xmax)] for i in (1 .. ymax)]
		for i in (1 .. ymax):
			T[i-1][0] = i
			for j in (xmin..xmax):
				T[i-1][j-xmin+1] = self.DeltaSet(j,i)
		T = [[''] + [j for j in (xmin..xmax)]] + T
		return T
	
	def HTMLDeltaSetTable(self, xmin, xmax, ymax):
		return table(self.DeltaSetTable(xmin,xmax,ymax))
	
	def Lambda(self):
		# max([i for i in self.irrs.values() if i < oo] + [i for i in self.heights.values() if i < oo])
		if self.__lambda == 0:
			self.__lambda = 1 + max([max([min(self.LengthSet(i,j)) for j in (1 .. min(self.ymax,self.LastReducibleAt(i)))]) for i in (1 .. self.frob+2*self.x0) if self.HeightOfRow(i) > 0])
		
		return self.__lambda
	
	# assumes factorizations are ordered (perhaps lexicographically)
	def __FactorDistance(self, fact1, fact2):
		if max(len(fact1),len(fact2)) - min(len(fact1),len(fact2)) > self.Lambda():
			return 0
		
		commons = 0
		i = 0
		j = 0
		
		while i < len(fact1) and j < len(fact2):
			if fact1[i] == fact2[j]:
				commons = commons + 1
				i = i + 1
				j = j + 1
			# comparison - key function?
			elif fact1[i] < fact2[j]:
				i = i + 1
			else:
				j = j + 1
		
		return max(len(fact1) - commons, len(fact2) - commons)
	
	def CatenaryDegree(self, x, n):
		#if (x,n) in self.catenaries.keys():
		#	return self.catenaries[(x,n)]
		
		if self.IsIrreducible(x,n):
			return 0
		
		# build the adjacency matrix
		facts = self.Factorizations(x,n)
    
		A = [[0 for f2 in facts] for f1 in facts]
		for i in range(0,len(facts)):
			for j in range(i+1,len(facts)):
				dist = self.__FactorDistance(facts[i],facts[j])
				# dists = dists.union(Set([dist]))
				A[i][j] = dist
				A[j][i] = dist
		
		dists = Set([d for l in A for d in l])
		dists = sorted([d for d in dists])
		
		# print dists
		
		# build the graph
		G = Graph(Matrix(A))
		
		# G.show()
		
		# prune the graph
		catenary = 0
		for catenary in dists:
			# if catenary == 0:
			# 	break
			
			if catenary == 0:
				continue
			
			Anew = [[(d if d <= catenary else 0) for d in l] for l in A]
			# print Anew
			'''
			for i in range(0, len(A)):
				for j in range(i, len(A)):
					if A[i][j] >= catenary:
						A[i][j] = 0
						A[j][i] = 0
			'''
			G = Graph(Matrix(Anew))
			# G.show()
			
			if G.is_connected():
				break
		
		# self.catenaries[(x,n)] = catenary
		return catenary
	
	def CatenaryDegreeTable(self, xmin, xmax, ymax):
		T = [[0] + [0 for j in (xmin..xmax)] for i in (1 .. ymax)]
		for i in (1 .. ymax):
			T[i-1][0] = i
			for j in (xmin..xmax):
				T[i-1][j-xmin+1] = self.CatenaryDegree(j,i)
		T = [[''] + [j for j in (xmin..xmax)]] + T
		return T
	
	def HTMLCatenaryDegreeTable(self, xmin, xmax, ymax):
		return table(self.CatenaryDegreeTable(xmin,xmax,ymax))
	
	def GammaPrime(self):
		primegens = []
		
		for i in (1 .. (self.frob+1)*2):
			if self.HeightOfRow(i) > 0 or self.HeightOfRow(i - self.s) > 0:
				primegens.append(i)
		
		primeret = NumericalSemigroup(primegens)
		return primeret
	
	def Plot(self, xmax = 0, ymax = 0):
		self.FindIrreducibles()
		
		if xmax == 0:
			xmax = self.xmax + 5
		if ymax == 0:
			ymax = self.ymax + 4
		
		p_irrs = []
		p_other = []
		for i in (1 .. xmax):
			if self.HeightOfRow(i) == 0:
				continue
			
			cutoff = min(self.FirstIrreducibleAt(i)-1, self.HeightOfRow(i), ymax)
			top = min(self.HeightOfRow(i), ymax)
			
			p_irrs.append((i,1))
			
			for j in (2 .. cutoff):
				p_other.append((i,j))
			
			for j in (cutoff+1 .. top):
				p_irrs.append((i,j))
		
		#P = point(points)
		P = point2d((0,0), figsize=(2*(xmax/ymax),4)) + point2d(p_irrs, rgbcolor='red', size=25) + point2d(p_other, rgbcolor='black',size=10) #, aspect_ratio=ymaximum/xmaximum)
		P.axes_range(0,xmax,0,ymax)
		return P
	
	# private helper functions
	def FindHeights(self):
		if self.heights != {}:
			return
		
		xmax = self.frob + 1
		
		self.heights = {i:0 for i in (1 .. xmax)}
		
		for i in self.heights.keys():
			if i not in self.gamma:
				continue
			
			j = i + self.s
			while j <= self.frob:
				if j not in self.gamma:
					break
				
				j = j + self.s
				self.heights[i] = self.heights[i] + 1
			
			if j > self.frob:
				self.heights[i] = oo
			
			if self.heights[i] > 0 and self.x0 == 0:
				self.x0 = i
			
			if self.heights[i] == oo and self.xf == 0:
				self.xf = i
			
		self.heights = {i:self.heights[i] for i in self.heights.keys() if i <= self.frob + self.x0}
	
	
	def FindIrreducibles(self):
		if self.irrs != {}:
			return
		
		self.FindHeights()
		
		xmax = self.frob + self.x0 + 1
		self.irrs = {i:2 for i in (1 .. xmax)}
		for i in (1 .. xmax):
			if self.HeightOfRow(i) == 0:
				self.irrs[i] = oo
				continue
			
			for j in (1 .. (i/2)+1):
				row1 = self.HeightOfRow(j)
				row2 = self.HeightOfRow(i-j)
				if row1 == 0 or row2 == 0:
					continue
				
				if row1 == oo or row2 == oo:
					self.irrs[i] = oo
					break
				
				self.irrs[i] = max(self.irrs[i], row1 + row2 + 1)
		
		self.xmax = self.frob + self.x0
		self.ymax = max([i for i in self.irrs.values() if i < oo] + [i for i in self.heights.values() if i < oo]) + 2
			
			
	    	







