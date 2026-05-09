using CertifiedHomotopyTracking

@variables x y
CC = AcbField(256)
G = [x^2 - 1, y^2 - 1]

F_com = [(1+1im)*x^2 + (3im)*y - 4, y^2 + 3]
H = straight_line_homotopy(F_com, G, [x, y]; CCRing=CC)
res = track_path(H, [CC(1), CC(-1)])
sol = solution(res)
evaluate_H(H, certified_region(res), CC(1))

F_com = [(1+1im)*x^2 + (3im)*y - 4im, y^2 + 3]
H = straight_line_homotopy(F_com, G, [x, y]; CCRing=CC)
res = track_path(H, [CC(1), CC(-1)])
sol = solution(res)
evaluate_H(H, certified_region(res), CC(1))

F_com = [(1+1im)*x^3/y + (3im)*y - 4, y^2 + 3]
H = straight_line_homotopy(F_com, G, [x, y]; CCRing=CC)
res = track_path(H, [CC(1), CC(-1)])
sol = certified_region(res)
evaluate_H(H, certified_region(res), CC(1))
solution(res)


using CertifiedHomotopyTracking

@variables x y 
PREC_BITS = 256
CC = AcbField(PREC_BITS)

F = [x^2 + 3*y - 4, y^2 + 3]
G = [x^2 - 1, y^2 - 1]
H = straight_line_homotopy(F, G, [x, y]; CCRing=CC)
start_point = [CC(1), CC(-1)]

evaluate_H(H, start_point, CC(0))


F_com = [(1+1im)*x^2 + 3*y - 4, y^2 + 3]
G = [x^2 - 1, y^2 - 1]
H = straight_line_homotopy(F_com, G, [x, y]; CCRing=CC)
start_point = [CC(1), CC(-1)]



evaluate_Jac(H, start_point, CC(0))
res = track_path(H, start_point)
sol = certified_region(res)
evaluate_H(H, certified_region(res), CC(1))