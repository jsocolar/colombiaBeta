db <- read.csv("/Users/jacobsocolar/Dropbox/CO_DBdata/DB_all_datasets.csv", sep = ";")
db2 <- db[, c("point", "trap", "day", "scientificName", "abundance")]
db2$point_day <- paste(db2$point, db2$day, sep = "__")
db2 <- db2[!is.na(db2$abundance), ] # trampas perdidas
db2 <- db2[db2$day != "", ]

db3 <- db2
db3$spd <- paste(db2$scientificName, db2$point_day, sep = "__")
db4 <- db3[duplicated(db3$spd) | duplicated(db3$spd, fromLast = T), ]

View(db4)

db2 <- db2[rownames(db2) != "209",] # remove extra uroxys coarctatus problem
db2 <- db2[rownames(db2) != "3357",]
db2 <- db2[rownames(db2) != "27640",]
db2 <- db2[rownames(db2) != "25096",]
db2 <- db2[rownames(db2) != "29548",]
db2 <- db2[rownames(db2) != "21280",]
db2 <- db2[rownames(db2) != "30185",]
db2 <- db2[rownames(db2) != "25097",]
db2 <- db2[rownames(db2) != "29549",]
db2 <- db2[rownames(db2) != "21281",]
db2 <- db2[rownames(db2) != "27642",]
db2 <- db2[rownames(db2) != "25098",]
db2 <- db2[rownames(db2) != "28914",]
db2 <- db2[rownames(db2) != "29550",]
db2 <- db2[rownames(db2) != "21282",]
db2 <- db2[rownames(db2) != "27643",]
db2 <- db2[rownames(db2) != "25099",]
db2 <- db2[rownames(db2) != "29551",]
db2 <- db2[rownames(db2) != "21283",]
db2 <- db2[rownames(db2) != "27644",]
db2 <- db2[rownames(db2) != "25100",]
db2 <- db2[rownames(db2) != "29552",]
db2 <- db2[rownames(db2) != "21284",]
db2 <- db2[rownames(db2) != "27645",]
db2 <- db2[rownames(db2) != "24465",]
db2 <- db2[rownames(db2) != "29553",]
db2 <- db2[rownames(db2) != "21285",]
db2 <- db2[rownames(db2) != "27648",]
db2 <- db2[rownames(db2) != "25104",]
db2 <- db2[rownames(db2) != "29556",]
db2 <- db2[rownames(db2) != "21288",]
db2 <- db2[rownames(db2) != "27649",]
db2 <- db2[rownames(db2) != "25105",]
db2 <- db2[rownames(db2) != "29557",]
db2 <- db2[rownames(db2) != "21289",]

# species <- unique(db2$scientificName)
# species <- species[species != ""]
# point_days <- unique(db2$point_day)

for (i in 2317:length(point_days)) { #seq_along(point_days)) {
  print(i)
  for (j in seq_along(species)) {
    L <- sum(db2$point_day == point_days[i] & db2$scientificName == species[j])
    if (L > 1) {stop(paste0("species ", j, "; point_day ", i))}
    if (L == 0) {
      db2 <- rbind(db2, data.frame(point = unique(db2$point[db2$point_day == point_days[i]]),
                                   trap = unique(db2$trap[db2$point_day == point_days[i]]),
                                   day = unique(db2$day[db2$point_day == point_days[i]]),
                                   scientificName = species[j],
                                   abundance = 0,
                                   point_day = point_days[i]))
    }
  }
}

db2[db2$scientificName == species[67] & db2$point_day == point_days[2317], ]
db2[db2$scientificName == species[140] & db2$point_day == point_days[2317], ]
db2[db2$scientificName == species[188] & db2$point_day == point_days[2317], ]
db2[db2$scientificName == species[220] & db2$point_day == point_days[2317], ]



db2 <- db2[db2$scientificName != "", ]
