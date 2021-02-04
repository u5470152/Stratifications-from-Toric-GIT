from ast import literal_eval


def ample_cone():
    M = matrix(fan.dim(), len(ray_list))
    relationList = []
    T = ToricVariety(fan)
    for cone in fan:
        for i in cone.ambient_ray_indices():
            M[:, i] = fan.rays()[i]
        for i in range(len(ray_list)):  # solve relation#
            if i not in cone.ambient_ray_indices():
                sol = M.solve_right(fan.rays()[i])
                relation = -sol
                relation[i] = 1
                if relation in relationList:  # avoid duplicates#
                    continue
                relationList.append(relation)
        M = matrix(fan.dim(), len(ray_list))
    for i in range(len(relationList)):
        relationList[i] = (0,) + tuple(relationList[i])  # to satisfy polyhedron format in sageMath#
    return Polyhedron(ieqs=relationList)



def reduce_ampleCone():
    l = []
    for normalVector in ample_cone_normal_vectors:
        reducedNormalVector = []
        for i in range(len(ray_list)):
            if i in choice:
                continue
            reducedNormalVector.append(normalVector[i])
        l.append((0,)+tuple(reducedNormalVector))
    return Polyhedron(ieqs = l)


def potential_stratum():  #computes the L-index of potential stratumans = []
	total_index = Set(range(len(ray_list)))
	for i in range(len(fan.primitive_collections())):
		for j in powerset(total_index - fan.primitive_collections()[i]):
			if list(j) in potential_stratum_index:
				continue
			potential_stratum_index .append(list(j))


def potential_worst_one_ps():
	tmp_matrix = R
	tmp_list = []
	vector_list=[]
	for i in range(len(potential_stratum_index)):
		for index in potential_stratum_index[i]:
			tmp_matrix = tmp_matrix.augment(standard_basis_matrix[:, index])
		if tmp_matrix.rank() == len(ray_list):
			tmp_matrix = R
			continue
		complement = tmp_matrix.left_kernel().basis_matrix()
		projection = matrix(complement.rank(), len(ray_list) - complement.rank()).augment(
			identity_matrix(complement.rank()))
		A = block_matrix([[tmp_matrix.transpose().row_space().basis_matrix()], [complement]])
		tmp_list.append(
			[complement.transpose() * projection * A.transpose().inverse()* coordinate_ample_cone_generic
				, potential_stratum_index[i]])
		tmp_matrix = R
	#group the worst one PS
	for x in tmp_list:
		x[0].set_immutable()
	values = set(x[0] for x in tmp_list)
	group_list = [[y, [x[1] for x in tmp_list if x[0] == y]] for y in values]
	for item in group_list:
		item = item.append((coordinate_ample_cone_generic.transpose() * item[0])[(0, 0)])
	return group_list


def is_proper_subspace(list1, list2):  #determine if the subsapce by list2 is a proper subspace defined by list1
	#list 2 should be longer
	m1 = R
	m2 = R
	for i in range(len(list1)):
		m1 = m1.augment(standard_basis_matrix[:,list1[i]])
	for i in range(len(list2)):
		m2 = m2.augment(standard_basis_matrix[:,list2[i]])
	space1 = m1.left_kernel().basis_matrix()
	space2 = m2.left_kernel().basis_matrix()
	return space1.nrows() > space2.nrows()


def is_co_dim_one_subspace(list_1, list_2): #second list produces a co-dim 1 subpace of list_1
	m1 = R
	m2 = R
	for i in range(len(list_1)):
		m1 = m1.augment(standard_basis_matrix[:, list_1[i]])
	for i in range(len(list_2)):
		m2 = m2.augment(standard_basis_matrix[:, list_2[i]])
	space1 = span(m1.left_kernel().basis_matrix(),QQ)
	space2 = span(m2.left_kernel().basis_matrix(),QQ)
	return(dim(space2) == dim(space1) - 1 and space2.is_subspace(space1))


def is_empty_type_one_wall(f): #takes in a linear expression
	coe_list =[]  #normal vector of f
	minus_coe_list = []
	for i in range(len(coordinate_ample_cone_reduced)):
		coe_list.append(f.coefficient({coordinate_ample_cone_reduced[i]:1}))
		minus_coe_list.append(-f.coefficient({coordinate_ample_cone_reduced[i]:1}))
	return dual_ample_cone_reduced.contains(coe_list) or dual_ample_cone_reduced.contains(minus_coe_list)



def type_one_wall():
	pre_list = []
	for s in potential_worst_one_ps_list:
		for t in potential_worst_one_ps_list:
			if is_co_dim_one_subspace(s[1][0],t[1][0]):
				m = R
				cutting_index = 0
				#pick out cutting index
				for i in range(len(s[1][0])):
					m = m.augment(standard_basis_matrix[:, s[1][0][i]])
				d = dim(span(m.left_kernel().basis_matrix(), QQ))
				for i in range(len(t[1][0])):
					m = m.augment(standard_basis_matrix[:, t[1][0][i]])
					if dim(span(m.left_kernel().basis_matrix(), QQ)) == d:
						continue
					cutting_index = t[1][0][i]
					break
				if is_empty_type_one_wall(s[0][(cutting_index,0)]):
					continue
				pre_list.append([s[0][(cutting_index, 0)],[s[1], t[1]]])
	#group type one walls
	values = set(x[0] for x in pre_list)
	return [[x, [y[1] for y in pre_list if y[0] == x]] for x in values]



def is_empty_type_two_wall(f, remaining_index, cross_polytope, cross_polytope_hrep):
	constraint_list = []
	#detect free variables
	non_free_variables = []
	for var in list(f.variables()):
		non_free_variables.append(ample_cone_coordinate_ring.gens().index(var))

	free_variables_indices = list(set(remaining_index).difference(set(non_free_variables)))
	if len(free_variables_indices) > 0:
		vertex_list = cross_polytope.vertices()
		vertices = list(list(vertex_list[i].vector()) for i in range(len(vertex_list)))
		#remember the index of free variables in the list of remaining_index
		remember = list(remaining_index.index(i) for i in free_variables_indices)
		#project cross_polytope by projecting vertices
		projected_vertices = []
		for item in vertices:
			projected_vertex = []
			for i in range(len(item)):
				if i in remember:
					continue
				projected_vertex.append(item[i])
			projected_vertices.append(projected_vertex)
		projected_cross_polytope = Polyhedron(vertices=projected_vertices)
		#get hrepresentation of projected_polytope

		constrained_variable = matrix(
			ample_cone_coordinate_ring.gens()[i] for i in list(set(remaining_index).difference(set(free_variables_indices))))\
			.transpose()
		for item in projected_cross_polytope.Hrepresentation():
			constraint_list.append(SR((matrix(item.A())*constrained_variable)[(0,0)]+item.b()))
		for i in free_variables_indices:
			remaining_index = remaining_index[:remaining_index.index(i)] + remaining_index[remaining_index.index(i)+1:]
		a = list(minimize_constrained(SR(f), constraint_list,
									list(0 for i in range(len(remaining_index)))))
		b = list(minimize_constrained(-SR(f), constraint_list,
									list(0 for i in range(len(remaining_index)))))

	else:
		for item in cross_polytope_hrep:
			constraint_list.append(SR(item))
		a = list(minimize_constrained(SR(f), constraint_list, list(0 for i in range(len(coordinate_ample_cone_reduced)-1))))
		b = list(minimize_constrained(-SR(f), constraint_list, list(0 for i in range(len(coordinate_ample_cone_reduced)-1))))
	a_point = list(0 for i in range(len(ray_list)))
	b_point = list(0 for i in range(len(ray_list)))
	counter = 0
	for j in range(len(remaining_index)):
		a_point[remaining_index[j]] = a[counter]
		b_point[remaining_index[j]] = b[counter]
		counter = counter + 1
	return f(a_point)*f(b_point) >= 0


def type_two_wall():
	pairs = []
	pre_list = []
	reduced_pre_list = []
	for s in potential_worst_one_ps_list:
		for t in potential_worst_one_ps_list:
			if ([t, s] not in pairs and is_proper_subspace(s[1][0], list(set(s[1][0]+t[1][0])))
			and is_proper_subspace(t[1][0], list(set(s[1][0]+t[1][0])))):
				pairs.append([s, t])
	for item in pairs:
		pre_list.append([item[0][2]-item[1][2], [item[0][1], item[1][1]]])

	#compute cross section
	A = matrix(coordinate_ample_cone_reduced).transpose()
	l = dual_ample_cone_reduced.rays()
	v = l[0]  # get something in the interior of the dual ample cone
	for i in range(len(l)):
		if i == 0:
			continue
		v = v + l[i]
	# cross section of the ample cone
	h = (matrix(v) * A)[(0, 0)]
	print('Cross Section: ')
	cross_section = h - (matrix(v) * matrix(ample_cone_reduced.rays()[0]).transpose())[(0, 0)]
	print(cross_section)
	var_index = Integer(raw_input('enter index of the variable to be substituted.'))
	expression_r = (-1 / cross_section.coefficient({ample_cone_coordinate_ring.gens()[var_index]: 1})) * \
				   (cross_section - cross_section.coefficient({ample_cone_coordinate_ring.gens()[var_index]: 1}) *
					ample_cone_coordinate_ring.gens()[var_index])
	#compute remaining index on the polytope
	index = set(choice).union({var_index})  # index of variables to be gone
	remaining_index = list(set(i for i in range(len(ray_list))).difference(index))
	#compute cross polytope
	cross_polytope_hrep = []
	for item in reduced_ample_cone_normal_vectors:
		cross_polytope_hrep.append((matrix(item) * A)[(0, 0)])
	for i in range(len(cross_polytope_hrep)):
		cross_polytope_hrep[i] = cross_polytope_hrep[i].substitute({ample_cone_coordinate_ring.gens()[var_index]: expression_r})

	cross_polytope = Polyhedron(ieqs=list((exp.constant_coefficient(),)
										  + tuple(
		exp.monomial_coefficient(ample_cone_coordinate_ring.gens()[i]) for i in remaining_index) for exp in
										  cross_polytope_hrep))

	for item in pre_list:
		reduced_eq = item[0].substitute({list(ample_cone_coordinate_ring.gens())[i]: 0 for i in list(choice)})
		if is_empty_type_two_wall(reduced_eq.substitute({ample_cone_coordinate_ring.gens()[var_index]: expression_r})
				, remaining_index, cross_polytope, cross_polytope_hrep):
			continue
		reduced_pre_list.append([reduced_eq, item[1]])
	plot_list = []
	plot_list.append(cross_polytope.plot())
	#show(cross_polytope.plot())
	for item in reduced_pre_list:
		wall_projection = SR(item[0].substitute({ample_cone_coordinate_ring.gens()[var_index]: expression_r}))
		plot_list.append(implicit_plot(wall_projection,(a2,0,1),(a3,0,2), region = lambda a2,a3: a2-a3+1>=0))
		show(cross_polytope.plot() + implicit_plot(wall_projection,(a2,0,1),(a3,0,2), region = lambda a2,a3: a2-a3+1>=0) )
	for item in type_one_wall_list:
		wall_projection = item[0].substitute({ample_cone_coordinate_ring.gens()[i]: 0 for i in choice})
		wall_projection = SR(wall_projection.substitute({ample_cone_coordinate_ring.gens()[var_index]: expression_r}))
		plot_list.append(implicit_plot(wall_projection,(a2,0,1),(a3,0,2), region = lambda a2,a3: a2-a3+1>=0, color = 'red'))
	show(sum(plot_list))


#computes the worst one ps for a single cone
def minimize_cone(l_index, character):
	"""
	:type character: vector
	"""
	power_set = list(powerset(l_index))
	tmp_list = []
	for item1 in power_set:
		for item2 in potential_worst_one_ps_list:
			for i in range(len(item2[1])):
				if item1 == item2[1][i]:
					if item2[0] not in tmp_list:
						tmp_list.append([item2[0], item2[2]])
	for item in tmp_list:
		item[0] = -item[0].substitute({ample_cone_coordinate_ring.gens()[i]:character[i] for i in range(len(character))})
		item[1] = -item[1].substitute({ample_cone_coordinate_ring.gens()[i]:character[i] for i in range(len(character))})
	#delete the one PS that are not in the cone
	unwanted = []
	for i in range(len(tmp_list)):
		for j in l_index:
			if tmp_list[i][0][j,:] < 0 and tmp_list[i] not in unwanted:
				unwanted.append(tmp_list[i])
	tmp_list = list(x for x in tmp_list if x not in unwanted)

	#insertion sort
	for i in range(1,len(tmp_list)):
		key = tmp_list[i]
		j = i-1
		while j >= 0 and key[1] < tmp_list[j][1]:
			tmp_list[j+1] = tmp_list[j]
			j = j-1
		tmp_list[j+1] = key
	tmp_list[0][0].set_immutable()
	return tuple(tmp_list[0])


#list out stratification
def test_stratification(character):
	tmp_list = []
	for item in potential_stratum_index:
		tmp_list.append((minimize_cone(item, character), item))
	#group the tori by worst one PS
	values = set(x[0] for x in tmp_list)
	return tuple((x, tuple(tuple(y[1]) for y in tmp_list if y[0] == x)) for x in values)


def is_ample(character):
	M = matrix([list(x) for x in ample_cone_normal_vectors])
	vector = matrix(len(ray_list),1)
	for i in range(len(ray_list)):
		vector[i,:] = character[i]
	result = M*vector
	for i in range(len(ample_cone_normal_vectors)):
		if result[(i, 0)] <= 0:
			return False
	return True

#the Stratification object is initiated by a character
class Stratification(object):
	def __init__(self, character):
		if not is_ample(character):
			raise ValueError('Not an ample divisor. Need help on figuring out an ample divisor? Type '
							 'ample_cone_normal_vectors')
		else:
			self.character = character
			self.strata = test_stratification(character)

	def is_worse(self, l1, l2):
		index1 = 0
		index2 = 0
		while self.strata[index1][1] != l1:
			index1 = index1 + 1
		while self.strata[index2][1] != l2:
			index2 = index2 + 1
		return self.strata[index1][0][1] > self.strata[index2][0][1]

	def ordering(self):
		return Poset((tuple(x[1] for x in self.strata),self.is_worse))

	def is_equivalent(self, strat):
		return self.ordering == strat.ordering()


potential_stratum_index = []  #L-index
potential_worst_one_ps_list = []
type_one_wall_list = []
type_two_wall_list = []
ample_cone_normal_vectors = []
reduced_ample_cone_normal_vectors = []
input_listRays = literal_eval(raw_input('Enter tuples for the rays, separated by commas.'))
input_listCones = literal_eval(raw_input('Enter the ambient indices for each cone.'))
ray_list = list(input_listRays)
coneList = list(input_listCones)
fan = Fan(cones=coneList, rays=ray_list)
standard_basis_matrix = identity_matrix(len(ray_list))
ample_cone_coordinate_ring = PolynomialRing(QQ,len(ray_list),"a")
ample_cone_coordinate_ring.inject_variables()# generic coordinate for ample cone#
coordinate_ample_cone_generic = matrix(ample_cone_coordinate_ring, len(ray_list),1)  #symbolic matrix
for i in range(len(ray_list)):
	coordinate_ample_cone_generic[i,:] = ample_cone_coordinate_ring.gens()[i]
coordinate_list = list(var('x_%d' % i) for i in range(len(ray_list)))
R = matrix(len(ray_list), fan.dim())  #the matrix M->Z^\{sigma}
for i in range(len(ray_list)):  # the matrix M->Z^{\Sigma}
	R[i, :] = fan.rays()[i]
#set up - ample cone, potential worst one parameter subgroups


potential_stratum() # enumerate all L-indicies
ample_cone_generic = ample_cone()
for item in ample_cone_generic.Hrepresentation():
	ample_cone_normal_vectors.append(item.A())
potential_worst_one_ps_list = potential_worst_one_ps()
choice = raw_input('What to do now? ample cone. stratification.')
if (choice == 'ample cone'):
	for item in ample_cone_normal_vectors:
		print(item)
if (choice == 'stratification'):
	for item in potential_worst_one_ps_list:
		print(item)
		print('\n')
choice = literal_eval(raw_input('choose a cone so that the basis of the divisor class group has basis indexed by'
								' the rays not in that cone'))  # returned choice is a tuple, representing the ambient indicies
# of the cone#
coordinate_ample_cone_reduced = list(ample_cone_coordinate_ring.gens()[i] for i in
								set(i for i in range(len(fan.rays()))).difference(set(choice)))
print('Reduced coordinate:')
print(coordinate_ample_cone_reduced)
ample_cone_reduced = reduce_ampleCone()
for item in ample_cone_reduced.Hrepresentation():
	reduced_ample_cone_normal_vectors.append(item.A())

print('The normal vectors of ample cone are: ')
print(reduced_ample_cone_normal_vectors)
dual_ample_cone_reduced = Cone(reduced_ample_cone_normal_vectors)
type_one_wall_list = type_one_wall()

#type_two_wall()
#testing data:(1,0),(0,1),(0,-1),(-1,1),(-1,2) \ (0,1),(0,2),(2,3),(3,4),(4,1)
#testing data:(3,4),(2,-3),(1,-2),(0,-1),(-1,1)\ (0,1),(0,2),(2,3),(3,4),(4,1)
#testing_data:(1,0,0),(0,1,0),(0,0,1),(0,-1,0),(1,1,1),(-1,-1,-1) \
# (0,1,4),(0,2,4),(1,2,4),(0,1,5),(0,3,5),(2,3,5),(1,2,5)