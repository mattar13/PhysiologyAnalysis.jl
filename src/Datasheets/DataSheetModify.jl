function expand_dates(datasheet)
    datasheet |> 
        @mutate(
            Date = Date(datasheet[1, "date_created"]), 
            Time = Time(datasheet[1, "date_created"])) |> 
    DataFrame
end