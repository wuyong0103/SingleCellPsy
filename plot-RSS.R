library(SCENIC)
rss <- readRDS("./Ex-L56NP/Ex-L56NP_rss.rds")
head(rss)
rss_filtered <- rss[apply(rss, 1, max) > 0.2, c("secondtrim", "thirdtrim", "infant", "juvenile", "youth", "midlife", "elder")]
pdf('RSS_plot_Ex-L56NP.pdf',width = 4, height = 8)
rssPlot <- plotRSS(rss_filtered, col.low = '#330066', col.mid = '#66CC66', col.high= '#FFCC33')
rssPlot$plot
dev.off()

rss <- readRDS("./In-SST/In-SST_rss.rds")
head(rss)
rss_filtered <- rss[apply(rss, 1, max) > 0.17, c("secondtrim", "thirdtrim", "infant", "kid", "child", "juvenile", "youth", "midlife", "elder")]
pdf('RSS_plot_In-SST.pdf',width = 4.78, height = 8)
rssPlot <- plotRSS(rss_filtered, col.low = '#330066', col.mid = '#66CC66', col.high= '#FFCC33')
rssPlot$plot
dev.off()
