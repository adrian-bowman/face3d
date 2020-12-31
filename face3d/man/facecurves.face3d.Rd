\name{facecurves.face3d}

\alias{facecurves.face3d}

\title{Find a planar, optimized, or smoothed path between anatomicaly pre-specified set of landmarks pairs along the surface of a face3d object}

\description{
This function finds a planar, optimized, or smoothed path between anatomicaly pre-specified set of landmarks pairs through the mesh of a face3d object.
}

\usage{
facecurves.face3d(face, curves, graphics = FALSE, monitor = FALSE, pcrv.path = FALSE,
                  start.pc = NA) 
}

\arguments{
\item{face}{A \code{face3d} object.}
\item{curves}{The specified curve to find.}
\item{graphics}{a logical value which controls whether diagnostic graphics are created as the function executes.}
\item{monitor}{a logical value which controls whether progress is monitored as the function executes.}
\item{pcrv.path}{a logical value which controls whether a principal curve or a planepath is used as a reference path.}
\item{start.pc}{an starting path used to initialise the principal curve algorithm.}
}

\details{
The 23 landmarks are as follows: pronasale (pn), alare crest L/R (acL, acR), subnasale (sn), sellion (se), nasion (n), exocanthion L/R (exL, exR), endocanthion L/R (enL, enR), tragion L/R (tL, tR), labiale superius (ls), crista philtri L/R (cphL, cphR), cheilion L/R (chL, chR), stomion (st), labiale inferius (li), sublabiale (sl), gnathion (gn), otobasion inferius L/R (oiL, oiR).

A planar, optimized, or smoothed path (curve) between anatomicaly pre-specified set of 24 landmarks pairs as follows: nasal root L/R (se-enL, enR-se), nasal boundary L/R (acL-enL, enR-acR), nasal bridge L/R (pn-acL, acR-pn), mid-line philtral (sn-ls), nasal base L/R (sn-acL, acR-sn), mid-line lip (chR-chL), lower lip L/R (li-chL, chR-li), upper lip L/R (cphL-chL, chR-cphR), philtrum ridge L/R (sn-cphL, cphR-sn), philtrum lip L/R (ls-cphL, cphR-ls), mid-line nasal profile (se-pn), mid-line nasal root (n-se), mid-line columella (pn-sn), mid-line upper lip (ls-st), mid-line bottom lip (st-li), mid-line mentolabial (li-sl), mid-line chin (sl-gn), nasolabial L/R (acL-chL, chR-acR), cheek-lip L/R (chL-oiL, oiR-chR), cheek-nose L/R (acL-tL, tR-acR), cheek-eye L/R (exL-tL, tR-exR), ear L/R (tL-oiL, tR-oiR), mandible L/R (gn-oiL, oiR-gn), Mandible (oiR-oiL), upper eye socket L/R (enL-exL,exR-enR), lower eye socket L/R (enL-exL,exR-enR), brow ridge L/R (n-tL, tR-n), Brow Ridge (tR-tL), open mouth curves (OM before lip curves).
}

\value{
Returning \code{shape} as \code{face3d} object whith a set of curves as a list concatenated to a matrix \code{shape$curves}.
}

\examples{
\dontrun{
   data(face)
   face <- facecurves.face3d(face, "nasal root right")
}}
