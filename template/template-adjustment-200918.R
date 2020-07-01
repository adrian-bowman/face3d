# Correct the triples on the lhs of the templates and add normals

load("~/owncloud/Shared/visualisation-shape/examples/data/template_male_symmetric.dmp")
load("~/owncloud/Shared/visualisation-shape/examples/data/template_female_symmetric.dmp")
template <- template.male
template <- template.female

trpls <- matrix(template$triples, ncol = 3, byrow = TRUE)
ind <- apply(trpls, 1, function(x) all(substr(rownames(template$coords)[x], 1, 1) == "L"))
trpls[ind, ] <- t(apply(trpls[ind, ], 1, rev))
template$triples <- c(t(trpls))
template <- normals.face3d(template)

# Find the normals for the template curves - an approximation for the moment
cnormals <- matrix(nrow = nrow(template$curves), ncol = 3)
for (i in 1:nrow(template$curves)) {
   ind <- which.min(edist.face3d(template$coords, template$curves[i, ]))
   cnormals[i, ] <- template$normals[ind, ]
}
template$cnormals <- cnormals

plot(template, display = "mesh", new = FALSE)
plot(template, display = "normals", add = TRUE)
spheres3d(template$curves, col = "red", radius = 0.2)
lns <- matrix(c(t(cbind(template$curves, template$curves + 5 * cnormals))), ncol = 3, byrow = TRUE)
segments3d(lns, col = "red")

template.male <- template
save(template.male, file = "~/owncloud/Face3D_0.1-1/Face3D/data/templateMale.rda")
template.female <- template
save(template.female, file = "~/owncloud/Face3D_0.1-1/Face3D/data/templateFemale.rda")
