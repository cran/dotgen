sim <- function()
{
    ## fetch data (already aligned by ID)
    res <- readRDS("rs208294_res.rds") # result
    gno <- readRDS("rs208294_gno.rds") # genotype
    cvr <- readRDS("rs208294_cvr.rds") # covariates

    ## condition on covariates
    cov.gno <- scp(cov(cbind(gno, cvr)), colnames(cvr))
    cor.gno <- cov2cor(cov.gno)
    stt.gno <- with(res, zvl(P, BETA)) # variant statistics

    ## -------- Combine statistics;  TQ approach --------
    ref <- tsq(stt.gno, cor.gno)
    ref <- do.call(.d, ref[c('mtd', 'ssq', 'pvl')])
    ## sum of quadratics = 50.53402; p-value = 0.001778937

    ## -------- Combine statistics; DOT approach --------
    alt <- dot(stt.gno, cor.gno)
    alt <- do.call(.d, alt[c('mtd', 'ssq', 'pvl')])
    ## sum of quadratics = 35.76306; p-value = 0.001132132

    rbind(ref, alt)
}
